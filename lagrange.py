# Copyright (C) 2016 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>

import jacobi
from matplotlib.backends.backend_pdf import PdfPages

if __name__=='__main__':
    with PdfPages('lagrange.pdf') as pdf_pages:
        fig=1
        for (Cj,limit) in zip([3.805,3.552,3.197,2.84011],[2,2,2,1]):    
            jacobi.plot_jacobi(fig=fig,Cj=Cj,pdf_pages=pdf_pages,limit=limit)
            fig+=1
