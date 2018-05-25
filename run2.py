from freqbands import *


if __name__ == '__main__':
    data_collect = glob.glob('*.rea')
    result = table_prod_4total(data_collect)
    print(result)
    table = np.sum(result,axis=1)/result.shape[1]
    print(table.shape)
    end_array =np.sum(result,axis=0)/len(data_collect)
    print(end_array.shape)
    end_array = np.reshape(end_array,(8,5))

    total_hypothesis_test = kruskal(table[:,0],table[:,1],table[:,2],table[:,3],
                        table[:,4,])
    print('total hypothesis: ',total_hypothesis_test)

    data_total = [table[:,0],table[:,1],table[:,2],table[:,3],
                        table[:,4,]]
    plt.figure()
    plt.boxplot(data_total)
    plt.title('Boxplots of total signal power')
    plt.show()
