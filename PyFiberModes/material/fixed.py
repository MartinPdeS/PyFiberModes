# This file is part of FiberModes.
#
# FiberModes is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FiberModes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FiberModes.  If not, see <http://www.gnu.org/licenses/>.

"""Module for fixed index material."""

from .material import Material


class Fixed(Material):

    """Fixed index material class.

    A material with a fixed index always have the same refractive index,
    whatever the wavelength is.

    """

    name = "Fixed index"
    info = "Fixed index"
    nparams = 1

    @classmethod
    def n(cls, wl, n):
        return n
