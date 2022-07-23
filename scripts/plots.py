

    plt.figure()
    for x,i in intra_subs_dis.items():
        if i>10:
            ind = d_complete.divergence==x
            plt.plot(sorted(d_complete.loc[ind,"numdate"]), np.linspace(0,1,i), label=f"n={x}")
    plt.legend()

    plt.figure()
    for x,i in intra_aaSubs_dis.items():
        if i>10:
            ind = d_complete.aaDivergence==x
            plt.plot(sorted(d_complete.loc[ind,"numdate"]), np.linspace(0,1,i), label=f"n={x}")
    plt.legend()

    plt.figure()
    for x,i in intra_subs_dis.items():
        if i>10:
            ind = d_complete.divergence==x
            plt.plot(sorted(d_complete.loc[ind,"numdate"]), np.arange(i), label=f"n={x}")
    plt.ylim(1,1000)
    plt.yscale('log')
    plt.legend()

    plt.figure()
    for x,i in intra_aaSubs_dis.items():
        if i>10:
            ind = d_complete.aaDivergence==x
            plt.plot(sorted(d_complete.loc[ind,"numdate"]), np.arange(i), label=f"n={x}")
    plt.ylim(1,1000)
    plt.yscale('log')
    plt.legend()

    plt.figure()
    for x,i in intra_geno.items():
        if i>100:
            ind = d_complete["intra_substitutions_str"]==x
            print(sorted(d_complete.loc[ind,"numdate"])[:3])
            plt.plot(sorted(d_complete.loc[ind,"numdate"]), np.arange(i)+1, label=f"n={x}")
    plt.ylim(0.5,1000)
    plt.yscale('log')
    plt.legend()


    plt.figure()
    plt.scatter(d.numdate, d.divergence)
