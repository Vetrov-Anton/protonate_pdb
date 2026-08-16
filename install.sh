#!/usr/bin/env bash
# Install protprep. Touches nothing outside its own virtual environment.
#
#   ./install.sh            - venv next to the sources (.venv) + protonate command
#   ./install.sh --pipx     - isolated install via pipx (command on PATH)
#   ./install.sh --user     - pip install --user (command in ~/.local/bin)
#
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mode="${1:---venv}"

case "$mode" in
  --pipx)
    command -v pipx >/dev/null || { echo "pipx not found: python3 -m pip install --user pipx"; exit 1; }
    pipx install --force "$here"
    echo
    echo "Done. Command: protonate --help"
    ;;
  --user)
    python3 -m pip install --user "$here"
    echo
    echo "Done. Command: protonate --help (needs ~/.local/bin on PATH)"
    ;;
  --venv|"")
    python3 -m venv "$here/.venv"
    "$here/.venv/bin/pip" install --quiet --upgrade pip
    # editable: source edits are picked up by both entry points
    "$here/.venv/bin/pip" install --quiet -e "$here"
    echo
    echo "Done. Command: $here/protonate --help"
    echo "(or source $here/.venv/bin/activate && protonate --help)"
    ;;
  *)
    echo "Unknown mode: $mode"; exit 2 ;;
esac
