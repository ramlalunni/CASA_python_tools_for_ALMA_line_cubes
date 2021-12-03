# combine_plots_in_latex.py

########################################################################
# Copyright (C) 2021  Ramlal Unnikrishnan

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# contact: ramlal.unnikrishnan@chalmers.se
########################################################################

import os
import subprocess

# function to combine all plots to single pdf file (by running 'combine_plots.tex' file present in pwd)
###----------------------------------------------------------------------------###
def combine_plots_in_latex(source):
    """combines plots produced by the functions for several line cubes (of the same source) into a single multi-page pdf file,
       with each page corresponding to the plots produced from a particualr line cube, by running the LaTeX script 'combine_plots.tex'.
       Requires pdflatex and pythontex to be installed.
       Requires 'combine_plots.tex' to exist in $PWD.
       Inputs: 1. source (str): name of the source corresponding to the input line cubes, Eg: 'IRAS_07454-7112'"""

    # removing files from previous compilation, if any
    os.system('rm -rf pythontex-files-combine_plots/ combine_plots.aux combine_plots.pytxcode combine_plots.log')

    # latex runs
    pdf_latex_run_1 = subprocess.run('pdflatex combine_plots.tex', shell=True, check=True, capture_output=False)
    if pdf_latex_run_1.returncode == 0:
        python_tex_run = subprocess.run('pythontex combine_plots.tex', shell=True, check=True, capture_output=False)
    if python_tex_run.returncode == 0:
        pdf_latex_run_2 = subprocess.run('pdflatex combine_plots.tex', shell=True, check=True, capture_output=False)
    if pdf_latex_run_2.returncode == 0:
        os.system('mv combine_plots.pdf %s_line_plots.pdf' %source)
        print('LaTeX compilation and run successful. All plots are saved to %s_line_plots.pdf' %source)

    os.system('rm -rf pythontex-files-combine_plots/ combine_plots.aux combine_plots.pytxcode combine_plots.log')
###----------------------------------------------------------------------------###
