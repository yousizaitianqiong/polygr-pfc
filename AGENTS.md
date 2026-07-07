# Agent Notes

## Windows Search Fallback

In this workspace on Windows, `rg.exe` may fail with `Access is denied`.
When that happens, fall back to native PowerShell search commands such as:

- `Get-ChildItem -Recurse -File | Select-String -Pattern ...`
- `Get-ChildItem -Recurse -Filter ...`

If the workspace contains both `polygr-pfc-fork` and `polygr-pfc-master`,
prefer `polygr-pfc-fork` when it has the active `.git` directory and recent
build artifacts, unless the user explicitly asks otherwise.
