# Примеры

* `fragment.pdb` — цепь A остатков 30–70 из 6CFO плюс кристаллические воды
  (~15 КБ). Хватает, чтобы за секунду проверить, что всё работает.
* `example_spec.txt` — файл задания для полного 6CFO (гомотетрамер): фиксация
  боковых цепей в обеих альфа-субъединицах, кэпы на цепях A и C.

Быстрая проверка (силовое поле берётся из установленного GROMACS по имени):

```bash
protonate -f examples/fragment.pdb -o /tmp/prot_demo --ph 7.0 \
    --ff amber99sb-ildn --fix A:31:p --fix A:63:HIE
```

Ожидаемо: Asp31 станет `ASH`, His63 — `HIE`, остальные группы — по PROPKA,
концы цепи заряжены. С `--run-pdb2gmx` тут же соберётся и топология.

Полный 6CFO (файл скачивается с RCSB, в репозиторий не входит):

```bash
wget https://files.rcsb.org/download/6CFO.pdb
protonate -f 6CFO.pdb -o prepared --spec examples/example_spec.txt \
    --ff amber99sb-ildn --drop-het --run-pdb2gmx
```
