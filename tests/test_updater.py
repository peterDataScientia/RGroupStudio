"""End-to-end updater contract tests without making network requests."""
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import updater

class _Response(io.BytesIO):
    def geturl(self):
        return "https://objects.githubusercontent.com/validated-installer.exe"

class UpdaterTests(unittest.TestCase):
    def _release(self, version="0.4.1", size=16):
        name=updater.SETUP_TEMPLATE.format(version=version)
        base=f"https://github.com/{updater.OWNER}/{updater.REPO}/releases/download/v{version}/"
        return {"tag_name":"v"+version,"draft":False,"prerelease":False,
                "assets":[
                  {"name":name,"browser_download_url":base+name,"size":size},
                  {"name":name+".sha256","browser_download_url":base+name+".sha256","size":91}]}

    def test_version_comparison(self):
        self.assertGreater(updater.version_parts("v0.2.10"),updater.version_parts("0.2.9"))
        with self.assertRaises(updater.UpdateError):
            updater.version_parts("0.2.9-malicious")

    def test_latest_release_found(self):
        payload=json.dumps(self._release()).encode()
        with patch("updater.read_url",return_value=payload):
            info=updater.latest_update()
        self.assertEqual(info["version"],"0.4.1")
        self.assertTrue(info["url"].endswith(info["filename"]))
        self.assertTrue(info["checksum_url"].endswith(".sha256"))

    def test_current_version_has_no_update(self):
        payload=json.dumps(self._release("0.4.0")).encode()
        with patch("updater.read_url",return_value=payload):
            self.assertIsNone(updater.latest_update())

    def test_reject_wrong_asset_source(self):
        release=self._release()
        release["assets"][0]["browser_download_url"]="https://example.org/malware.exe"
        with patch("updater.read_url",return_value=json.dumps(release).encode()):
            with self.assertRaises(updater.UpdateError):
                updater.latest_update()

    def test_reject_missing_hash(self):
        release=self._release()
        release["assets"].pop()
        with patch("updater.read_url",return_value=json.dumps(release).encode()):
            with self.assertRaises(updater.UpdateError):
                updater.latest_update()

    def test_download_verifies_checksum_and_bytes(self):
        payload=b"verified installer"
        release=self._release(size=len(payload))
        with patch("updater.read_url",return_value=json.dumps(release).encode()):
            info=updater.latest_update()
        digest=hashlib.sha256(payload).hexdigest()
        manifest=(digest+"  "+info["filename"]+"\n").encode()
        with patch("updater.read_url",return_value=manifest), \
             patch("updater.urllib.request.urlopen",return_value=_Response(payload)):
            file=updater.download_installer(info)
        self.assertEqual(file.read_bytes(),payload)
        file.unlink()
        file.parent.rmdir()

    def test_reject_checksum_mismatch(self):
        payload=b"damaged installer"
        release=self._release(size=len(payload))
        with patch("updater.read_url",return_value=json.dumps(release).encode()):
            info=updater.latest_update()
        manifest=("0"*64+"  "+info["filename"]+"\n").encode()
        with patch("updater.read_url",return_value=manifest), \
             patch("updater.urllib.request.urlopen",return_value=_Response(payload)):
            with self.assertRaises(updater.UpdateError):
                updater.download_installer(info)

    def test_auto_check_setting(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"update_settings.json"
            with patch("updater.settings_path",return_value=p):
                self.assertTrue(updater.get_auto_check())
                updater.set_auto_check(False)
                self.assertFalse(updater.get_auto_check())
                updater.set_auto_check(True)
                self.assertTrue(updater.get_auto_check())

if __name__=="__main__":
    unittest.main()
