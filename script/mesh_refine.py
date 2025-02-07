
def refine_by_splitting(x,y,z):
    nx1 = len(x[0][0])
    ny1 = len(x[0])
    nz1 = len(x)

    x2 = []; y2 = []; z2 = []
    for k in range(nz1):
        X2D = []; Y2D = []; Z2D = []
        for j in range(ny1-1):
            X1D = []; Y1D = []; Z1D = []
            for i in range(nx1-1):
                X1D.append(x[k][j][i])
                X1D.append(0.5*(x[k][j][i]+x[k][j][i+1]))
                Y1D.append(y[k][j][i])
                Y1D.append(0.5*(y[k][j][i]+y[k][j][i+1]))
                Z1D.append(z[k][j][i])
                Z1D.append(0.5*(z[k][j][i]+z[k][j][i+1]))
            X1D.append(x[k][j][i+1])
            Y1D.append(y[k][j][i+1])
            Z1D.append(z[k][j][i+1])
            X2D.append(X1D)
            Y2D.append(Y1D)
            Z2D.append(Z1D)
            X1D = []; Y1D = []; Z1D = []
            for i in range(nx1-1):
                X1D.append(0.5*(x[k][j][i]+x[k][j+1][i]))
                X1D.append(0.5*(x[k][j][i]+x[k][j+1][i])+0.5*(-x[k][j][i]+x[k][j][i+1]))
                Y1D.append(0.5*(y[k][j][i]+y[k][j+1][i]))
                Y1D.append(0.5*(y[k][j][i]+y[k][j+1][i])+0.5*(-y[k][j][i]+y[k][j][i+1]))
                Z1D.append(0.5*(z[k][j][i]+z[k][j+1][i]))
                Z1D.append(0.5*(z[k][j][i]+z[k][j+1][i])+0.5*(-z[k][j][i]+z[k][j][i+1]))
            X1D.append(0.5*(x[k][j][i+1]+x[k][j+1][i+1]))
            Y1D.append(0.5*(y[k][j][i+1]+y[k][j+1][i+1]))
            Z1D.append(0.5*(z[k][j][i+1]+z[k][j+1][i+1]))
            X2D.append(X1D)
            Y2D.append(Y1D)
            Z2D.append(Z1D)
        X1D = []; Y1D = []; Z1D = []
        for i in range(nx1-1):
            X1D.append(0.5*(x[k][j+1][i]+x[k][j+1][i]))
            X1D.append(0.5*(x[k][j+1][i]+x[k][j+1][i])+0.5*(-x[k][j+1][i]+x[k][j+1][i+1]))
            Y1D.append(0.5*(y[k][j+1][i]+y[k][j+1][i]))
            Y1D.append(0.5*(y[k][j+1][i]+y[k][j+1][i])+0.5*(-y[k][j+1][i]+y[k][j+1][i+1]))
            Z1D.append(0.5*(z[k][j+1][i]+z[k][j+1][i]))
            Z1D.append(0.5*(z[k][j+1][i]+z[k][j+1][i])+0.5*(-z[k][j+1][i]+z[k][j+1][i+1]))
        X1D.append(0.5*(x[k][j+1][i+1]+x[k][j+1][i+1]))
        Y1D.append(0.5*(y[k][j+1][i+1]+y[k][j+1][i+1]))
        Z1D.append(0.5*(z[k][j+1][i+1]+z[k][j+1][i+1]))
        X2D.append(X1D)
        Y2D.append(Y1D)
        Z2D.append(Z1D)
        x2.append(X2D)
        y2.append(Y2D)
        z2.append(Z2D)

    return x2,y2,z2