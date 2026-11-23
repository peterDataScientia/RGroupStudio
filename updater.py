"""Safe, installer-based RGroupStudio update checks using official GitHub Releases.

Updates never replace a running executable and are only installed with user approval.
The installer is downloaded to a temporary folder and SHA-256 checked before launch.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import urllib.request
from urllib.parse import urlparse

APP_VERSION = "0.3.0"
OWNER = "peterDataScientia"
REPO = "RGroupStudio"
API_URL = f"https://api.github.com/repos/{OWNER}/{REPO}/releases/latest"
SETUP_TEMPLATE = "RGroupStudio_Setup_v{version}_Windows_x64.exe"
MAX_SIZE = 250 * 1024 * 1024
USER_AGENT = f"RGroupStudio-Updater/{APP_VERSION}"

class UpdateError(Exception):
    pass

def version_parts(value):
    version = str(value).strip().lstrip("vV")
    if not re.fullmatch(r"[0-9]+[.][0-9]+[.][0-9]+", version):
        raise UpdateError("Unsupported version: " + str(value))
    return tuple(int(x) for x in version.split("."))

def settings_path():
    base = Path(os.getenv("LOCALAPPDATA") or Path.home()) / "RGroupStudio"
    return base / "update_settings.json"

def get_auto_check():
    try:
        content = json.loads(settings_path().read_text(encoding="utf-8"))
        return content.get("auto_check", True) is not False
    except (OSError, ValueError, TypeError):
        return True

def set_auto_check(value):
    path = settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"auto_check": bool(value)}, indent=2), encoding="utf-8")

def read_url(url, limit=1048576, timeout=12):
    req = urllib.request.Request(url, headers={
        "User-Agent": USER_AGENT, "Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        content = response.read(limit + 1)
    if len(content) > limit:
        raise UpdateError("Update response is too large")
    return content

def _trusted_asset(asset, filename, version):
    if not isinstance(asset, dict) or asset.get("name") != filename:
        raise UpdateError("Missing official release asset: " + filename)
    url = asset.get("browser_download_url", "")
    expected = f"https://github.com/{OWNER}/{REPO}/releases/download/v{version}/{filename}"
    if url != expected:
        raise UpdateError("Unexpected release download URL")
    return url

def latest_update():
    """Return release info only when newer than installed version, else None."""
    try:
        release = json.loads(read_url(API_URL).decode("utf-8"))
    except (ValueError, UnicodeError) as exc:
        raise UpdateError("Invalid GitHub release response") from exc
    if not isinstance(release, dict) or release.get("draft") or release.get("prerelease"):
        raise UpdateError("No stable GitHub release available")
    tag = release.get("tag_name", "")
    version = str(tag).lstrip("vV")
    if version_parts(version) <= version_parts(APP_VERSION):
        return None
    filename = SETUP_TEMPLATE.format(version=version)
    checksum_name = filename + ".sha256"
    assets = {a.get("name"): a for a in release.get("assets", []) if isinstance(a, dict)}
    installer = assets.get(filename)
    checksum = assets.get(checksum_name)
    url = _trusted_asset(installer, filename, version)
    checksum_url = _trusted_asset(checksum, checksum_name, version)
    size = installer.get("size")
    if not isinstance(size, int) or size <= 0 or size > MAX_SIZE:
        raise UpdateError("Invalid installer size")
    return {
        "version": version, "filename": filename, "url": url,
        "checksum_url": checksum_url, "size": size,
        "release_url": f"https://github.com/{OWNER}/{REPO}/releases/tag/v{version}",
    }

def expected_sha256(info):
    blob = read_url(info["checksum_url"], limit=512).decode("ascii").strip()
    # Format: 64 hex digits followed by whitespace and exact filename.
    parts = blob.split()
    if len(parts) != 2 or not re.fullmatch("[a-fA-F0-9]{64}", parts[0]):
        raise UpdateError("Malformed checksum manifest")
    if parts[1].lstrip("*") != info["filename"]:
        raise UpdateError("Checksum filename mismatch")
    return parts[0].lower()

def download_installer(info, progress=None):
    expected = expected_sha256(info)
    folder = Path(tempfile.mkdtemp(prefix="RGroupStudio-update-"))
    dest = folder / info["filename"]
    temp_file = folder / (info["filename"] + ".part")
    total = 0
    digest = hashlib.sha256()
    try:
        req = urllib.request.Request(info["url"], headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=45) as response, temp_file.open("wb") as handle:
            origin = response.geturl()
            parsed = urlparse(origin)
            if parsed.scheme != "https" or not parsed.hostname:
                raise UpdateError("Untrusted download redirect")
            while True:
                block = response.read(256 * 1024)
                if not block:
                    break
                total += len(block)
                if total > MAX_SIZE or total > info["size"]:
                    raise UpdateError("Installer exceeds declared size")
                digest.update(block)
                handle.write(block)
                if progress:
                    progress(total, info["size"])
        if total != info["size"]:
            raise UpdateError("Installer download incomplete")
        if digest.hexdigest().lower() != expected:
            raise UpdateError("Installer SHA-256 verification failed")
        temp_file.replace(dest)
        return dest
    except Exception:
        temp_file.unlink(missing_ok=True)
        dest.unlink(missing_ok=True)
        try:
            folder.rmdir()
        except OSError:
            pass
        raise
