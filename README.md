# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                               |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|--------------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| PyFiberModes/analysis.py                           |       70 |       10 |       34 |        1 |     83.65% |209-216, 232-233 |
| PyFiberModes/coordinates.py                        |       69 |        3 |        4 |        2 |     93.15% |48, 124, 145 |
| PyFiberModes/fiber.py                              |      276 |       42 |       52 |        8 |     81.71% |62, 113, 223, 249-252, 272, 274, 279, 309, 336, 416-420, 461, 530-544, 614-618, 789, 891-895, 910-918, 937-938, 955-956 |
| PyFiberModes/field.py                              |      250 |       85 |       50 |        7 |     63.33% |122, 238-240, 261-265, 317-318, 372-375, 462-471, 516-525, 564-566, 602-611, 651-654, 696-697, 739-740, 804-814, 844-851, 882-889, 920-927, 962-974, 1076-1080, 1140, 1157-1159 |
| PyFiberModes/fundamentals.py                       |       79 |       26 |       24 |        6 |     59.22% |36-38, 64-67, 148-153, 187, 257-260, 304, 314, 325-332 |
| PyFiberModes/loader.py                             |       52 |        5 |       16 |        3 |     88.24% |47, 82, 184-186 |
| PyFiberModes/materials.py                          |       55 |        8 |       12 |        2 |     82.09% |29-31, 41, 91, 99-101 |
| PyFiberModes/models.py                             |       23 |        3 |       10 |        3 |     81.82% |30, 34, 63 |
| PyFiberModes/optimization.py                       |       40 |        1 |        6 |        1 |     95.65% |       164 |
| PyFiberModes/plotting.py                           |       27 |        6 |       10 |        4 |     72.97% |46-47, 53, 67, 73, 75 |
| PyFiberModes/propagation.py                        |       54 |        2 |       12 |        2 |     93.94% |   37, 174 |
| PyFiberModes/services.py                           |       56 |        9 |       10 |        2 |     83.33% |34, 68-71, 75, 112-113, 157 |
| PyFiberModes/solver/base\_solver.py                |       86 |       20 |       28 |       10 |     71.93% |54, 80-82, 130-\>133, 136, 142-143, 148, 154-\>158, 160-166, 207, 212, 254-262, 302-304 |
| PyFiberModes/solver/multilayer/effective\_index.py |      222 |       75 |       66 |        9 |     66.67% |56, 62, 71, 74-81, 100-101, 131-\>134, 135, 156-157, 165, 198-229, 264-276, 355, 455, 478-535, 715, 734 |
| PyFiberModes/solver/results.py                     |       15 |        1 |        2 |        1 |     88.24% |        56 |
| PyFiberModes/solver/three\_layer/cutoff.py         |      214 |      190 |       72 |        0 |      8.39% |30-31, 63-82, 97-126, 149-188, 203-223, 240-259, 276-290, 307-323, 340-356, 373-384, 399-424, 440-471, 486-503 |
| PyFiberModes/solver/two\_layer/cutoff.py           |       54 |       19 |       14 |        3 |     64.71% |55-56, 86-102, 139-152, 172-173 |
| PyFiberModes/solver/two\_layer/effective\_index.py |      198 |       94 |       32 |        8 |     52.17% |58-61, 87, 115, 136, 138, 140, 141-\>exit, 218, 318-337, 375-402, 422-479, 497, 577-579, 602-610, 691-696 |
| PyFiberModes/source.py                             |       17 |        5 |        0 |        0 |     70.59% |35, 47, 59, 70, 81 |
| PyFiberModes/stepindex.py                          |      135 |        4 |       28 |        4 |     95.09% |43, 45, 47, 95 |
| **TOTAL**                                          | **2099** |  **608** |  **490** |   **76** | **68.10%** |           |

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