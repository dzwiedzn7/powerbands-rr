from freqbands import *
"""
one of two main runs (use case)

after running this run(xD) you will obtain:
-csv with data how diffrent windows function and frequencies affect PSD of frequency
band of rr intervals signal

-result of kruskal test (testing if diffrent windows really affects PSD of frequency bands)
-boxplots of given PSDs
"""

if __name__ == '__main__':
    data_collect = glob.glob('*.rea')
    result = table_prod_4freqrange(data_collect)
    print(result)
    table = np.sum(result,axis=1)/result.shape[1]
    print(table.shape)
    end_array =np.sum(result,axis=0)/len(data_collect)
    print(end_array.shape)
    end_array = np.reshape(end_array,(8,15))

    hypothesis_test = kruskal(table[:,0,0],table[:,1,0],table[:,2,0],table[:,3,0],
                             table[:,4,0],nan_policy='omit')
    print ('vlf_hypothesis_test: ',hypothesis_test)

    hypothesis_test = kruskal(table[:,0,1],table[:,1,1],table[:,2,1],table[:,3,1],
                                     table[:,4,1])
    print ('lf_hypothesis_test: ',hypothesis_test)

    hypothesis_test = kruskal(table[:,0,2],table[:,1,2],table[:,2,2],table[:,3,2],
                                    table[:,4,2],nan_policy='omit')
    print ('hf_hypothesis_test: ',hypothesis_test)


    data_vlf = [table[:,0,0],table[:,1,0],table[:,2,0],table[:,3,0],
                         table[:,4,0]]
    plt.figure()
    plt.boxplot(data_vlf)
    plt.title('Boxplots of vlf band power')
    plt.show()
    data_lf = [table[:,0,1],table[:,1,1],table[:,2,1],table[:,3,1],
                         table[:,4,1]]


    plt.figure()
    plt.boxplot(data_lf)
    plt.title('Boxplots of lf band power')
    plt.show()


    data_hf = [table[:,0,2],table[:,1,2],table[:,2,2],table[:,3,2],
                        table[:,4,2]]
    plt.figure()
    plt.boxplot(data_hf)
    plt.title('Boxplots of hf band power')
    plt.show()
