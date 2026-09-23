# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                               |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|--------------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyFiberModes/analysis.py                           |       67 |       10 |       32 |        1 |     82.83% |206-213, 229-230 |
| PyFiberModes/coordinates.py                        |       69 |        3 |        4 |        2 |     93.15% |48, 124, 145 |
| PyFiberModes/fiber.py                              |      266 |       41 |       46 |        7 |     81.41% |62, 113, 223, 249-252, 272, 274, 279, 309, 336, 416-420, 461, 530-544, 614-618, 881-885, 900-908, 927-928, 945-946 |
| PyFiberModes/field.py                              |      250 |       85 |       50 |        7 |     63.33% |122, 238-240, 261-265, 317-318, 372-375, 462-471, 516-525, 564-566, 602-611, 651-654, 696-697, 739-740, 804-814, 844-851, 882-889, 920-927, 962-974, 1076-1080, 1140, 1157-1159 |
| PyFiberModes/fundamentals.py                       |       75 |       27 |       20 |        4 |     56.84% |35-37, 63-66, 141-149, 183, 253-256, 300, 319-326 |
| PyFiberModes/loader.py                             |       52 |        5 |       16 |        3 |     88.24% |47, 82, 184-186 |
| PyFiberModes/materials.py                          |       55 |        8 |       12 |        2 |     82.09% |29-31, 41, 91, 99-101 |
| PyFiberModes/models.py                             |       23 |        3 |       10 |        3 |     81.82% |30, 34, 63 |
| PyFiberModes/optimization.py                       |       40 |        1 |        6 |        1 |     95.65% |       164 |
| PyFiberModes/plotting.py                           |       27 |        6 |       10 |        4 |     72.97% |46-47, 53, 67, 73, 75 |
| PyFiberModes/propagation.py                        |       54 |        2 |       12 |        2 |     93.94% |   37, 174 |
| PyFiberModes/services.py                           |       51 |        7 |        8 |        2 |     84.75% |29, 61-62, 66, 103-104, 148 |
| PyFiberModes/solver/base\_solver.py                |       86 |       20 |       28 |       10 |     71.93% |54, 80-82, 130-\>133, 136, 142-143, 148, 154-\>158, 160-166, 207, 212, 254-262, 302-304 |
| PyFiberModes/solver/multilayer/effective\_index.py |      222 |       75 |       66 |        9 |     66.67% |43, 49, 58, 61-68, 87-88, 118-\>121, 122, 143-144, 152, 185-216, 251-263, 342, 442, 465-522, 702, 721 |
| PyFiberModes/solver/results.py                     |       15 |        1 |        2 |        1 |     88.24% |        56 |
| PyFiberModes/solver/three\_layer/cutoff.py         |      214 |      190 |       72 |        0 |      8.39% |25-26, 53-72, 87-116, 139-178, 193-213, 230-249, 266-280, 297-313, 330-346, 363-374, 389-414, 430-461, 476-493 |
| PyFiberModes/solver/two\_layer/cutoff.py           |       54 |       19 |       14 |        3 |     64.71% |45-46, 76-92, 129-142, 162-163 |
| PyFiberModes/solver/two\_layer/effective\_index.py |      198 |       94 |       32 |        8 |     52.17% |45-48, 74, 102, 123, 125, 127, 128-\>exit, 205, 305-324, 362-389, 409-466, 484, 563-565, 585-593, 670-675 |
| PyFiberModes/source.py                             |       17 |        5 |        0 |        0 |     70.59% |35, 47, 59, 70, 81 |
| PyFiberModes/stepindex.py                          |      135 |        4 |       28 |        4 |     95.09% |43, 45, 47, 95 |
| **TOTAL**                                          | **2063** |  **606** |  **472** |   **73** | **67.61%** |           |

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