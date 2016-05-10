import matplotlib
import time
import warnings
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from functools import reduce
from matplotlib import figure
import numpy as np
import operator
from numpy import *
from numpy import linalg as LA
from numpy.linalg import inv, norm
from tkinter import *
from tkinter.ttk import *

warnings.filterwarnings("ignore")


# GaussSeidel Solution Method
def Gauss_Seidel(A, B, x, tol):
    # A, x, b - matrices defining the problem
    # tol - absolute tolerance for convergence. (Can be adjusted)
    Y = np.tril(A)
    Z = A - Y
    Max_Iter = 10000
    iters = 0
    for i in range(Max_Iter):
        iters += 1
        x_prev = x
        x = np.dot(np.linalg.inv(Y), B - np.dot(Z, x))
        error = norm(x - x_prev) / norm(x)
        if error < tol:
            return x, iters


# Creates Mesh That Defines System
def Vars(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint):
    meshx = np.linspace(xrange[0], xrange[1], n)
    meshy = np.linspace(yrange[0], yrange[1], n)
    S = np.zeros((n, n))
    if isPoint[0] == 1:  # Checks to See if Source is a Point/(s)
        for i in range(0, len(Source)):
            S[isPoint[i + 1][1]][isPoint[i + 1][0]] = Source[i]
    elif isPoint[0] == 0:  # Checks to See if Source is Uniform
        S = Source[0] * np.ones((n, n))
    D = 1 / (3 * Sig_T)
    Diff = D * np.ones((n + 1, n + 1))
    SigA = Sig_A * np.ones((n + 1, n + 1))
    delx = abs(meshx[1] - meshx[0]) * ones(n + 1)
    dely = abs(meshy[1] - meshy[0]) * ones(n + 1)
    return meshx, meshy, S, SigA, Diff, delx, dely


# Finite Volume Method Solution-------------------------------------------------------------------------------------
def FinVol(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol):
    start_time = time.time()
    # Getting Mesh and Variables
    meshx, meshy, Source, SigA, D, delx, dely = Vars(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint)

    # defining Matrix Coefficients
    ALeft = np.zeros((n, n))
    ARight = np.zeros((n, n))
    ATop = np.zeros((n, n))
    ABot = np.zeros((n, n))
    ACenter = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            ALeft[i][j] = -(D[i][j] * dely[j] + D[i][j + 1] * dely[j + 1]) / (2 * delx[i])
            ARight[i][j] = -(D[i][j] * dely[j] + D[i + 1][j + 1] * dely[j + 1]) / (2 * delx[i + 1])
            ABot[i][j] = -(D[i][j] * delx[i] + D[i + 1][j] * delx[i + 1]) / (2 * dely[j])
            ATop[i][j] = -(D[i][j + 1] * delx[i] + D[i + 1][j + 1] * delx[i + 1]) / (2 * dely[j + 1])
            ACenter[i][j] = SigA[i][j] - (ALeft[i][j] + ARight[i][j] + ABot[i][j] + ATop[i][j])

    # Boundary Conditions
    for j in range(1, int(n)):
        ARight[0][j] = 0
        ABot[0][j] = 0  # Vacuum On Left Side
        ATop[0][j] = 0
        ACenter[0][j] = 1
        Source[j][0] = 0
    for j in range(0, int(n) - 1):
        ALeft[j][0] = 0
        ARight[j][0] = 0  # Vacuum On Bottom
        ATop[j][0] = 0
        ACenter[j][0] = 1
        Source[0][j] = 0
    for j in range(1, int(n)):
        ALeft[n - 1][j] = -(D[n - 2][j] * dely[j] + D[n - 2][j - 1] * dely[j - 1]) / (2 * delx[n - 2])
        ABot[n - 1][j] = -D[n - 2][j] * delx[j] / (2 * dely[j])  # Reflecting on Right
        ATop[n - 1][j] = -D[n - 2][j - 1] * delx[j] / (2 * dely[j + 1])
        ACenter[n - 1][j] = SigA[n - 1][j] - (ALeft[n - 1][j] + ABot[n - 1][j] + ATop[n - 1][j])
    for j in range(0, int(n) - 1):
        ARight[j][n - 1] = -D[i][n - 1] * dely[n - 1] / (2 * delx[j])  # Reflecting on Top
        ABot[j][n - 1] = -(D[j][n - 1] * delx[j] + D[i + 1][n - 1] * delx[i + 1]) / (2 * dely[n - 1])
        ALeft[j][n - 1] = -D[j + 1][n - 1] * dely[n - 1] / (2 * delx[j + 1])
        ACenter[j][n - 1] = SigA[j][n - 1] - (ARight[j][n - 1] + ABot[j][n - 1] + ALeft[j][n - 1])
    Center = []

    # Assignment of Coefficients to nxnxn Matrix
    for i in range(n):
        Center.append([] * (n))
        Center[i] = zeros((n, n))
    for j in range(0, n):
        for i in range(1, n - 1):
            Center[j][i][i] = ACenter[i][j]
            Center[j][i][i + 1] = ARight[i][j]
            Center[j][i][i - 1] = ALeft[i][j]
        Center[j][0][0] = ACenter[0][j]
        Center[j][0][1] = ARight[0][j]
        Center[j][n - 1][n - 1] = ACenter[n - 1][j]
        Center[j][n - 1][n - 2] = ALeft[n - 1][j]
    Top = []
    for i in range(n):
        Top.append([] * n)
        Top[i] = zeros((n, n))
    for j in range(0, n):
        for i in range(0, n):
            Top[i][j][j] = ATop[j][i]
    Bot = []
    for i in range(n):
        Bot.append([] * n)
        Bot[i] = zeros((n, n))
    for j in range(n):
        for i in range(0, n):
            Bot[i][j][j] = ABot[j][i]

    # Pulling Coefficients to (n**2Xn**2) to Matrix
    A = np.zeros((n ** 2, n ** 2))
    c0 = []
    c1 = []
    c2 = []
    b = []
    t = []
    N = n ** 2
    for i in range(n):
        temp1 = (diag(Center[i]))
        for j in range(n):
            c1.append(temp1[j])
    for i in range(n):
        temp0 = (diag(Center[i], -1))
        temp2 = (diag(Center[i], 1))
        if i > 0:
            temp0 = insert(temp0, 0, 0)
            temp2 = insert(temp2, 0, 0)
        for j in range(len(temp0)):
            c0.append(temp0[j])
            c2.append(temp2[j])
    for i in range(1, n):
        tempb = (diag(Bot[i]))
        for j in range(len(tempb)):
            b.append(tempb[j])
    for i in range(0, n - 1):
        tempt = (diag(Top[i]))
        for j in range(len(tempt)):
            t.append(tempt[j])
    for i in range(len(c1)):
        A[i][i] = c1[i]
    for i in range(len(c0)):
        A[i + 1][i] = c0[i]
        A[i][i + 1] = c2[i]
    for i in range(len(b)):
        A[i + n][i] = b[i]
    for i in range(len(t)):
        A[i][i + n] = t[i]

    # Creating Source Matrix and Initial Guess Matrix (Our Guess Will Just Be Our Source)
    Source = np.reshape(Source, (N, 1))
    guess = Source

    # Performing a Gauss-Seidel Iteration Method
    Flux, iters = Gauss_Seidel(A, Source, guess, tol)
    # Print TOC and Iterations
    print("Time of Completion: %s seconds" % (time.time() - start_time))
    print("Number of Iterations:" + str(iters))
    return Flux, iters, meshx, meshy


# Plotting Functions-------------------------------------------------------------------------------------------------
# Creating a Surface Plot of Flux Solution
def Surf_Plot():
    n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol = input_vars()
    Flux, iters, meshx, meshy = (FinVol(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol))
    XX, YY = meshgrid(meshx, meshy)
    flux = np.reshape(Flux, (n, n))
    fig = plt.figure(figsize=(16, 9))
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(XX, YY, flux, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5, label="Flux[n/cm^2/s")
    ax.set_xlabel('X Position [cm]')
    ax.set_ylabel('Y Position [cm]')
    ax.set_zlabel('Flux [n/cm^2/s')
    plt.title('Surface Map of Flux Solutions to 2D Fixed Source Diffusion Equation', size=22)
    plt.show()


# Views from X-Z Axis and X-Z Axis
def Side_View():
    n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol = input_vars()
    Flux, iters, meshx, meshy = (FinVol(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol))
    xs = np.linspace(xrange[0], xrange[1], len(Flux))
    ys = np.linspace(yrange[0], yrange[1], len(Flux))
    fig = plt.figure(figsize=(16, 9))
    ax1 = plt.subplot(122)
    ax2 = plt.subplot(121)
    ax1.plot(xs, Flux, 'b.', label="X-Z View")
    ax2.plot(ys, Flux, 'r.', label="Y-Z View")
    ax1.legend()
    ax2.legend()
    ax1.set_xlabel('X Position [cm]')
    ax1.set_ylabel('Flux [n/cm^2/s')
    ax2.set_xlabel('Y Position [cm]')
    ax2.set_ylabel('Flux [n/cm^2/s')
    fig.suptitle("Side Views of Flux Solutions to 2D Fixed Source Diffusion Equation", fontsize="x-large")
    plt.show()


# Gives a Heat Map of the Flux
def Heat_Map():
    n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol = input_vars()
    Flux, iters, meshx, meshy = (FinVol(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol))
    XX, YY = meshgrid(meshx, meshy)
    flux = np.reshape(Flux, (n, n))
    fig = plt.figure(figsize=(16, 9))
    plt.pcolor(XX, YY, flux, cmap=cm.coolwarm)
    plt.xlabel('X Position [cm]')
    plt.ylabel('Y Position [cm]')
    plt.title("Heat Map of Flux Solutions to 2D Fixed Source Diffusion Equation", fontsize="x-large")
    plt.colorbar(label="Flux[n/cm^2/s")
    plt.show()
# End Plotting Functions--------------------------------------------------------------------------------------------


# Creating Gui-----------------------------------------------------------------------------------------------------

# Closes Gui
def close_window():
    master.destroy()


# Enables and Disable Source Location if Source is Uniform
def dis():
    if var.get() == 1:
        SourceL.config(state=DISABLED)
    else:
        SourceL.config(state=NORMAL)


# Inputs All Data Stored In Text Boxes
def input_vars():
    N = 0
    N = int(n.get())
    xrange = [-1, 1]
    yrange = [-1, 1]
    xrange[0] = float(xmin.get())
    xrange[1] = float(xmax.get())
    yrange[0] = float(ymin.get())
    yrange[1] = float(ymax.get())
    s = Source.get().split(',')
    S = []
    Sig_t = 0.45
    Sig_a = 0.10
    for i in s:
        S.append(int(i))
    Sig_t = float(Sig_T.get())
    Sig_a = float(Sig_A.get())
    SL = []
    sl = SourceL.get().split(',')
    for i in sl:
        SL.append(int(i))
    isPoint = [0]
    SL = reshape(SL, (len(SL) / 2, 2))
    if var.get() == 0:  # Checks to See if Source is a Point/(s)
        isPoint[0] = 1
        for i in range(0, len(SL)):
            isPoint.append([SL[i][0], SL[i][1]])
    else:
        isPoint[0] = 0
    TOL = float(tol.get())
    return N, xrange, yrange, S, Sig_t, Sig_a, isPoint, TOL


# Export all our data to a TXT file
def Export():
    n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol = input_vars()
    Flux, iters, meshx, meshy = (FinVol(n, xrange, yrange, Source, Sig_T, Sig_A, isPoint, tol))
    np.savetxt(export_name.get(), Flux, delimiter=',')


# tkinter portion
master = Tk()
style = Style()
style.configure("BW.TLabel", foreground="black", background="white")
style.theme_use('classic')
Tk.iconbitmap(master, default='icon.ico')
Tk.wm_title(master, "2D Diffusion Equation Solver                       By:Will Kable")
Label(master, text="Mesh Size").grid(row=1, column=3)
Label(master, text="Lower X Boundary [cm]").grid(row=2, column=3)
Label(master, text="Upper X Boundary [cm]").grid(row=3, column=3)
Label(master, text="Lower Y Boundary [cm]").grid(row=4, column=3)
Label(master, text="Upper Y Boundary [cm]").grid(row=5, column=3)
Label(master, text="Source(s) [#/cm^2/s]").grid(row=1, column=5)
Label(master, text="Source Mesh Location").grid(row=2, column=5)
Label(master, text="Sigma Transport [1/cm]").grid(row=3, column=5)
Label(master, text="Sigma Absorption [1/cm]").grid(row=4, column=5)
Label(master, text="Tolerance").grid(row=5, column=5)

n = Entry(master)
n.insert(0, 30)
n.grid(row=1, column=4)
xmin = Entry(master)
xmin.insert(0, -1)
xmin.grid(row=2, column=4)
xmax = Entry(master)
xmax.insert(0, 1)
xmax.grid(row=3, column=4)
ymin = Entry(master)
ymin.insert(0, -1)
ymin.grid(row=4, column=4)
ymax = Entry(master)
ymax.insert(0, 1)
ymax.grid(row=5, column=4)
Source = Entry(master)
Source.insert(0, '10000,1000')
Source.grid(row=1, column=6)
SourceL = Entry(master)
SourceL.insert(0, '10,20,2,2')
SourceL.grid(row=2, column=6)
Sig_T = Entry(master)
Sig_T.insert(0, 0.45)
Sig_T.grid(row=3, column=6)
Sig_A = Entry(master)
Sig_A.insert(0, 0.1)
Sig_A.grid(row=4, column=6)
tol = Entry(master)
tol.insert(0, 0.0000001)
tol.grid(row=5, column=6)
export_name = Entry(master)
export_name.insert(0, 'Out.txt')
export_name.grid(row=5, column=1)
var = IntVar()
Checkbutton(master, text="Is the Source Uniform", variable=var, command=dis).grid(row=6, column=6, sticky=W, pady=4)
Button(master, text='Close', command=close_window).grid(row=6, column=1, sticky=W, pady=4)
Button(master, text='Surface Plot', command=Surf_Plot).grid(row=1, column=1, sticky=W, pady=4)
Button(master, text='Side Plot', command=Side_View).grid(row=2, column=1, sticky=W, pady=4)
Button(master, text='Heat Map', command=Heat_Map).grid(row=3, column=1, sticky=W, pady=4)
Button(master, text='Export Data', command=Export).grid(row=4, column=1, sticky=W, pady=4)

master.geometry('%dx%d+%d+%d' % (680, 260, 0, 0))
master.mainloop()
# end GUI Creation-------------------------------------------------------
