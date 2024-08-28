# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                      |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------ | -------: | -------: | -------: | -------: | ------: | --------: |
| PyFiberModes/\_\_future\_\_.py            |       48 |       48 |        6 |        0 |      0% |     4-106 |
| PyFiberModes/\_\_init\_\_.py              |       10 |        2 |        0 |        0 |     80% |     11-12 |
| PyFiberModes/\_version.py                 |       11 |        2 |        2 |        1 |     77% |       5-6 |
| PyFiberModes/factory.py                   |       47 |        0 |       10 |        0 |    100% |           |
| PyFiberModes/fiber.py                     |      210 |       38 |       54 |       12 |     76% |55->54, 65->64, 75->74, 85->84, 95->94, 105, 117, 139-141, 197-200, 203->202, 228-230, 233->232, 250->249, 260-264, 285->284, 297, 300->299, 348-362, 418-426, 580->583, 593->592, 606-613, 628-637, 640-648 |
| PyFiberModes/field.py                     |      215 |      168 |       54 |        0 |     18% |43, 56, 59-61, 64-67, 81-97, 111-127, 141-153, 167-182, 196-212, 226-239, 253-262, 276-291, 305-322, 336-352, 366-376, 390-409, 423-439, 453-460, 474-487, 501-516, 525-533, 545-551, 563-575, 587-594, 603, 606-610, 613-620 |
| PyFiberModes/fundamentals.py              |       70 |       38 |       22 |        2 |     41% |23-25, 40-45, 76, 106-118, 181-184, 212-246 |
| PyFiberModes/mode.py                      |       28 |        9 |       10 |        0 |     55% |30-35, 41-45 |
| PyFiberModes/mode\_instances.py           |       21 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/solver/\_\_init\_\_.py       |        1 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/solver/base\_solver.py       |       70 |       14 |       26 |        9 |     74% |23, 38->41, 44, 50-51, 56, 62->66, 68-74, 113, 118, 166-167 |
| PyFiberModes/solver/mlsif/\_\_init\_\_.py |        1 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/solver/mlsif/neff.py         |      225 |      171 |       64 |        7 |     22% |16-17, 37-38, 41, 44, 50-53, 56-63, 75-82, 85-86, 102-103, 121-182, 198-229, 261-273, 293, 311-355, 373, 395-452, 491-511, 514-536, 539-578, 581 |
| PyFiberModes/solver/ssif/\_\_init\_\_.py  |        2 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/solver/ssif/cutoff.py        |       55 |       19 |       14 |        3 |     65% |33-34, 48-64, 75-88, 108-109 |
| PyFiberModes/solver/ssif/neff.py          |      194 |      103 |       32 |        8 |     46% |32-35, 58, 83, 101, 103, 105, 106->exit, 141, 184-198, 232-251, 286-313, 330-387, 404, 470-472, 489-497, 565-570 |
| PyFiberModes/solver/tlsif/\_\_init\_\_.py |        1 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/solver/tlsif/cutoff.py       |      215 |      189 |       74 |        0 |      9% |23-24, 29-48, 51-80, 98-133, 145-165, 179-198, 212-226, 240-256, 270-286, 300-311, 314-339, 343-374, 377-394 |
| PyFiberModes/stepindex.py                 |      150 |       96 |       28 |        2 |     36% |26, 50, 121-122, 173-240, 249-289, 299-335 |
| PyFiberModes/tools/\_\_init\_\_.py        |        0 |        0 |        0 |        0 |    100% |           |
| PyFiberModes/tools/directories.py         |       16 |       16 |        4 |        0 |      0% |      4-39 |
| PyFiberModes/tools/utils.py               |       11 |        0 |        4 |        0 |    100% |           |
| PyFiberModes/wavelength.py                |       49 |       17 |       30 |        8 |     53% |30, 32, 36, 39-52, 57->56, 64->63, 73->72, 80->79, 84 |
|                                 **TOTAL** | **1650** |  **930** |  **434** |   **52** | **40%** |           |


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