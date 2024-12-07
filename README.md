# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------ | -------: | -------: | -------: | -------: | ------: | --------: |
| PyFiberModes/\_\_future\_\_.py      |       48 |       48 |        6 |        0 |      0% |     4-106 |
| PyFiberModes/coordinates.py         |       67 |        7 |        2 |        0 |     87% |46, 88, 100, 117-120 |
| PyFiberModes/factory.py             |       47 |        0 |        8 |        0 |    100% |           |
| PyFiberModes/fiber.py               |      210 |       38 |       30 |        1 |     78% |105, 117, 139-141, 197-200, 228-230, 260-264, 297, 348-362, 418-426, 580->583, 606-613, 628-637, 640-648 |
| PyFiberModes/field.py               |      238 |      184 |       60 |        0 |     18% |24-30, 33, 75, 116-122, 142-144, 164-166, 212-229, 274-292, 330-343, 388-404, 449-465, 499-512, 546-555, 589-604, 639-656, 691-707, 737-747, 777-796, 826-842, 873-880, 911-924, 959-974, 995-1003, 1023-1029, 1051-1063, 1086-1090, 1096, 1118-1135 |
| PyFiberModes/fundamentals.py        |       67 |       37 |       20 |        1 |     40% |29-31, 57-60, 131-137, 203-206, 241-275 |
| PyFiberModes/loader.py              |       46 |        4 |       14 |        2 |     90% |42, 162-164 |
| PyFiberModes/mode.py                |       28 |        9 |        8 |        0 |     53% |30-35, 41-45 |
| PyFiberModes/mode\_instances.py     |       21 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/solver/base\_solver.py |       70 |       14 |       26 |        9 |     74% |23, 38->41, 44, 50-51, 56, 62->66, 68-74, 113, 118, 166-167 |
| PyFiberModes/solver/mlsif/neff.py   |      225 |      171 |       64 |        7 |     22% |16-17, 37-38, 41, 44, 50-53, 56-63, 75-82, 85-86, 102-103, 121-182, 198-229, 261-273, 293, 311-355, 373, 395-452, 491-511, 514-536, 539-578, 581 |
| PyFiberModes/solver/ssif/cutoff.py  |       55 |       19 |       14 |        3 |     65% |33-34, 48-64, 75-88, 108-109 |
| PyFiberModes/solver/ssif/neff.py    |      194 |      103 |       32 |        8 |     46% |32-35, 58, 83, 101, 103, 105, 106->exit, 141, 184-198, 232-251, 286-313, 330-387, 404, 470-472, 489-497, 565-570 |
| PyFiberModes/solver/tlsif/cutoff.py |      215 |      189 |       72 |        0 |      9% |23-24, 29-48, 51-80, 98-133, 145-165, 179-198, 212-226, 240-256, 270-286, 300-311, 314-339, 343-374, 377-394 |
| PyFiberModes/stepindex.py           |      127 |       80 |       22 |        1 |     35% |54, 112, 239-306, 334-372, 410-444 |
| PyFiberModes/tools/utils.py         |       11 |        0 |        4 |        0 |    100% |           |
| PyFiberModes/wavelength.py          |       39 |        1 |       12 |        1 |     96% |        53 |
|                           **TOTAL** | **1708** |  **904** |  **394** |   **33** | **43%** |           |


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