# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyFiberModes/analysis.py            |       64 |        9 |       30 |        1 |     82.98% |203-209, 219-220 |
| PyFiberModes/coordinates.py         |       69 |        3 |        4 |        2 |     93.15% |48, 124, 145 |
| PyFiberModes/fiber.py               |      265 |       41 |       46 |        7 |     81.35% |61, 112, 222, 248-251, 271, 273, 278, 308, 335, 415-419, 460, 529-543, 613-617, 880-884, 899-907, 926-927, 944-945 |
| PyFiberModes/field.py               |      269 |       88 |       58 |       10 |     64.53% |130, 246-248, 269-273, 325-326, 380-383, 470-479, 524-533, 572-574, 610-619, 659-662, 704-705, 747-748, 812-822, 852-859, 890-897, 928-935, 970-982, 1084-1088, 1148, 1165-1167, 1200, 1207, 1210 |
| PyFiberModes/fundamentals.py        |       68 |       26 |       20 |        3 |     55.68% |29-31, 57-60, 132-138, 204-207, 247, 264-271 |
| PyFiberModes/loader.py              |       46 |        4 |       14 |        2 |     90.00% |41, 161-163 |
| PyFiberModes/models.py              |       23 |        3 |       10 |        3 |     81.82% |30, 34, 63 |
| PyFiberModes/optimization.py        |       40 |        1 |        6 |        1 |     95.65% |       164 |
| PyFiberModes/propagation.py         |       54 |        2 |       12 |        2 |     93.94% |   37, 174 |
| PyFiberModes/services.py            |       54 |        7 |        8 |        2 |     85.48% |31, 69-70, 74, 111-112, 156 |
| PyFiberModes/solver/base\_solver.py |       86 |       18 |       28 |       10 |     73.68% |54, 80-82, 130-\>133, 136, 142-143, 148, 154-\>158, 160-166, 207, 212, 302-304 |
| PyFiberModes/solver/mlsif/neff.py   |      194 |       99 |       58 |        7 |     50.00% |42, 48, 57, 60-67, 99-\>102, 103-104, 120-121, 143-174, 209-221, 247, 268-312, 331, 354-411, 591, 610 |
| PyFiberModes/solver/results.py      |       15 |        1 |        2 |        1 |     88.24% |        56 |
| PyFiberModes/solver/ssif/cutoff.py  |       55 |       19 |       14 |        3 |     65.22% |50-51, 81-97, 134-147, 167-168 |
| PyFiberModes/solver/ssif/neff.py    |      194 |       94 |       32 |        8 |     51.33% |49-52, 78, 106, 127, 129, 131, 132-\>exit, 170, 267-286, 324-351, 371-428, 446, 525-527, 547-555, 632-637 |
| PyFiberModes/solver/tlsif/cutoff.py |      215 |      190 |       72 |        0 |      8.71% |30-31, 58-77, 92-121, 144-180, 195-215, 232-251, 268-282, 299-315, 332-348, 365-376, 391-416, 432-463, 478-495 |
| PyFiberModes/source.py              |       17 |        5 |        0 |        0 |     70.59% |35, 47, 59, 70, 81 |
| PyFiberModes/stepindex.py           |      135 |        4 |       28 |        4 |     95.09% |43, 45, 47, 95 |
| **TOTAL**                           | **1952** |  **614** |  **446** |   **66** | **65.47%** |           |

5 files skipped due to complete coverage.


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/MartinPdeS/PyFiberModes/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/MartinPdeS/PyFiberModes/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2FMartinPdeS%2FPyFiberModes%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.