# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyFiberModes/\_\_future\_\_.py      |       48 |       48 |        6 |        0 |      0.00% |     3-175 |
| PyFiberModes/analysis.py            |       64 |        9 |       30 |        1 |     82.98% |203-209, 219-220 |
| PyFiberModes/coordinates.py         |       67 |        3 |        4 |        2 |     92.96% |48, 123, 143 |
| PyFiberModes/fiber.py               |      216 |       43 |       30 |        1 |     76.42% |154, 169, 195-198, 264-267, 299-301, 327-331, 372, 446-460, 525-533, 739-\>742, 768-775, 793-802, 817-825, 844-845, 862-863 |
| PyFiberModes/field.py               |      283 |       93 |       64 |       10 |     65.13% |46-52, 62, 181, 297-299, 320-324, 376-377, 431-434, 521-530, 575-584, 623-625, 661-670, 710-713, 755-756, 798-799, 863-873, 903-910, 941-948, 979-986, 1021-1033, 1135-1139, 1199, 1216-1218, 1250-\>1252, 1258, 1261 |
| PyFiberModes/fundamentals.py        |       68 |       37 |       20 |        1 |     40.91% |29-31, 57-60, 132-138, 204-207, 242-276 |
| PyFiberModes/loader.py              |       46 |        4 |       14 |        2 |     90.00% |41, 161-163 |
| PyFiberModes/mode.py                |       25 |        9 |        8 |        0 |     48.48% |51-56, 76-80 |
| PyFiberModes/optimization.py        |       40 |        1 |        6 |        1 |     95.65% |       164 |
| PyFiberModes/propagation.py         |       54 |        2 |       12 |        2 |     93.94% |   37, 174 |
| PyFiberModes/solver/base\_solver.py |       70 |       14 |       26 |        9 |     73.96% |51, 95-\>98, 101, 107-108, 113, 119-\>123, 125-131, 173, 178, 231-232 |
| PyFiberModes/solver/mlsif/neff.py   |      225 |      126 |       64 |        7 |     44.98% |24-25, 58, 64, 73, 76-83, 115-\>118, 119-120, 136-137, 159-220, 240-271, 306-318, 344, 365-409, 428, 451-508, 688, 707 |
| PyFiberModes/solver/ssif/cutoff.py  |       55 |       19 |       14 |        3 |     65.22% |50-51, 81-97, 134-147, 167-168 |
| PyFiberModes/solver/ssif/neff.py    |      194 |      103 |       32 |        8 |     46.46% |49-52, 78, 106, 127, 129, 131, 132-\>exit, 170, 216-230, 267-286, 324-351, 371-428, 446, 525-527, 547-555, 632-637 |
| PyFiberModes/solver/tlsif/cutoff.py |      215 |      190 |       72 |        0 |      8.71% |30-31, 58-77, 92-121, 144-180, 195-215, 232-251, 268-282, 299-315, 332-348, 365-376, 391-416, 432-463, 478-495 |
| PyFiberModes/source.py              |       17 |        5 |        0 |        0 |     70.59% |35, 47, 59, 70, 81 |
| PyFiberModes/stepindex.py           |      124 |        2 |       22 |        1 |     97.95% |   53, 115 |
| **TOTAL**                           | **1871** |  **708** |  **428** |   **48** | **59.37%** |           |

2 files skipped due to complete coverage.


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