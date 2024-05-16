import dipoleq
import matplotlib.pyplot as plt
import numpy as np

m = dipoleq.MACHINE('../../Testing/beta1.in')
m.LHname=""
dipoleq.solve(m)

plt.contour(m.PsiGrid.R, m.PsiGrid.Z, m.PsiGrid.Psi)
plt.show()
