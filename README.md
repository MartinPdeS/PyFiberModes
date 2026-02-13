# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyFiberModes/\_\_future\_\_.py      |       48 |       48 |        6 |        0 |      0.00% |     4-106 |
| PyFiberModes/coordinates.py         |       67 |        3 |        4 |        2 |     92.96% |46, 120, 140 |
| PyFiberModes/fiber.py               |      209 |       38 |       30 |        1 |     77.82% |144, 156, 178-180, 239-242, 270-272, 302-306, 346, 417-431, 487-495, 691->694, 717-724, 739-748, 751-759 |
| PyFiberModes/field.py               |      232 |      119 |       60 |        6 |     45.55% |24-30, 33, 164-166, 226-227, 289-292, 388-404, 449-465, 506-508, 546-555, 597-600, 639-656, 691-707, 737-747, 777-796, 826-842, 873-880, 911-924, 959-974, 1086-1090, 1096, 1130, 1133 |
| PyFiberModes/fundamentals.py        |       68 |       37 |       20 |        1 |     40.91% |30-32, 58-61, 133-139, 205-208, 243-277 |
| PyFiberModes/loader.py              |       46 |        4 |       14 |        2 |     90.00% |42, 162-164 |
| PyFiberModes/mode.py                |       25 |        9 |        8 |        0 |     48.48% |30-35, 41-45 |
| PyFiberModes/solver/base\_solver.py |       70 |       14 |       26 |        9 |     73.96% |23, 38->41, 44, 50-51, 56, 62->66, 68-74, 113, 118, 166-167 |
| PyFiberModes/solver/mlsif/neff.py   |      225 |      126 |       64 |        7 |     44.98% |16-17, 38, 44, 53, 56-63, 81->84, 85-86, 102-103, 121-182, 198-229, 261-273, 293, 311-355, 373, 395-452, 576, 581 |
| PyFiberModes/solver/ssif/cutoff.py  |       55 |       19 |       14 |        3 |     65.22% |33-34, 48-64, 75-88, 108-109 |
| PyFiberModes/solver/ssif/neff.py    |      194 |      103 |       32 |        8 |     46.46% |32-35, 58, 83, 101, 103, 105, 106->exit, 141, 184-198, 232-251, 286-313, 330-387, 404, 470-472, 489-497, 565-570 |
| PyFiberModes/solver/tlsif/cutoff.py |      215 |      190 |       72 |        0 |      8.71% |22-23, 28-47, 50-79, 97-133, 145-165, 179-198, 212-226, 240-256, 270-286, 300-311, 314-339, 343-374, 377-394 |
| PyFiberModes/source.py              |       17 |        5 |        0 |        0 |     70.59% |21, 33, 45, 56, 67 |
| PyFiberModes/stepindex.py           |      124 |        2 |       22 |        1 |     97.95% |   54, 112 |
| PyFiberModes/tools/directories.py   |       16 |       16 |        4 |        0 |      0.00% |      4-39 |
| **TOTAL**                           | **1682** |  **733** |  **384** |   **40** | **53.29%** |           |

3 files skipped due to complete coverage.


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