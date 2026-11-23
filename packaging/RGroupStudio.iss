#define AppName "RGroupStudio"
#define AppVersion "0.3.0"
[Setup]
AppId={{70F2B1E1-63A9-41A0-A6BE-4D77824479E7}
AppName={#AppName}
AppVersion={#AppVersion}
SetupIconFile=..\assets\RGroupStudio.ico
DefaultDirName={localappdata}\Programs\RGroupStudio
DefaultGroupName=RGroupStudio
OutputDir=..\release
OutputBaseFilename=RGroupStudio_Setup_v{#AppVersion}_Windows_x64
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
PrivilegesRequired=lowest
UninstallDisplayIcon={app}\RGroupStudio.exe
DisableProgramGroupPage=yes
CloseApplications=yes
[Files]
Source: "..\dist\RGroupStudio\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
[Icons]
Name: "{autoprograms}\RGroupStudio"; Filename: "{app}\RGroupStudio.exe"
Name: "{autodesktop}\RGroupStudio"; Filename: "{app}\RGroupStudio.exe"; Tasks: desktopicon
[Tasks]
Name: "desktopicon"; Description: "Create Desktop Shortcut"; GroupDescription: "Additional shortcuts:"; Flags: checkedonce
[Run]
Filename: "{app}\RGroupStudio.exe"; Description: "Launch RGroupStudio"; Flags: nowait postinstall skipifsilent
