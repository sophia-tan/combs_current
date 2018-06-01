import matplotlib.pyplot as plt
    sortedactivation = sorted(activation)
        new_order = [np.where(activation==i)[0][0] for ix,i in enumerate(sortedactivation)]
    calc = [calc[ix] for ix in new_order]
    AAs = [AAs[ix] for ix in new_order]
    activation = [activation[ix] for ix in new_order]
    
    def best_fit_slope(xs,ys):
        xs,ys=np.array(xs),np.array(ys)
        m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys)) / \
                 ((np.mean(xs)*np.mean(xs)) - np.mean(xs*xs)))
        b = np.mean(ys) - m*np.mean(xs)
        print('m',m)
        print('b',b)
        return m,b
    
    m,b=best_fit_slope(calc,activation)
    regression_line = [(m*x)+b for x in calc]
    ##x=[0.001,11.875]
    ##y = [0.18405186840597873,0.67299933809206247]
    #x=[-1.1,2]
    #y=[0.17258020395350002, 0.582692403831]

    corr = stats.pearsonr(activation,calc)
    print(corr[0],corr[1])
    
    
    for label, x, y in zip(AAs, calc,activation):
        text(label,x,y)
    
    ax.scatter(calc,activation,marker='o',edgecolors='gray',facecolors=col,
            lw=1)

plt.tight_layout()
plt.show()
