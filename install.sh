#!/usr/bin/env bash
# Установка protprep. Ничего кроме собственного venv не трогает.
#
#   ./install.sh            - venv рядом с исходниками (.venv) + команда protonate
#   ./install.sh --pipx     - изолированная установка через pipx (команда в PATH)
#   ./install.sh --user     - pip install --user (команда в ~/.local/bin)
#
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mode="${1:---venv}"

case "$mode" in
  --pipx)
    command -v pipx >/dev/null || { echo "pipx не найден: python3 -m pip install --user pipx"; exit 1; }
    pipx install --force "$here"
    echo
    echo "Готово. Команда: protonate --help"
    ;;
  --user)
    python3 -m pip install --user "$here"
    echo
    echo "Готово. Команда: protonate --help (нужен ~/.local/bin в PATH)"
    ;;
  --venv|"")
    python3 -m venv "$here/.venv"
    "$here/.venv/bin/pip" install --quiet --upgrade pip
    # editable: правки в исходниках сразу видны обеим точкам входа
    "$here/.venv/bin/pip" install --quiet -e "$here"
    echo
    echo "Готово. Команда: $here/protonate --help"
    echo "(или source $here/.venv/bin/activate && protonate --help)"
    ;;
  *)
    echo "Неизвестный режим: $mode"; exit 2 ;;
esac
