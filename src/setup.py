
# Setting up global variables

#in_files_path = '/home/aazlm/Documents/Git/fits_files/input/'
#out_files_path = '/home/aazlm/Documents/Git/fits_files/output/'

in_files_path = '../input/'
out_files_path = '../output/'

frequencies = ['100', '143', '217', '353', '545', '857']
K_to_MJy_coef = [244.1, 371.74, 483.69, 287.45, 58.04, 2.27]

# ILC components
SZ_comp = [-9.2634830362821635e2, -9.3931697105364112e2, 4.5363163361568719e1,
           1.461967447538203e3, 9.2614619777207065e2, 1.3059846139207055e2]

dust_comp = [0.00544493, 0.01867971,  0.07611707,
             0.36487403,  1.32332848,  4.1680303]

cold_dust_comp = [0.05917835,   0.19296831,   0.71168128,
                  2.76059493,   7.08596294, 11.7202328 ]

hot_dust_comp = [2.83225192e-03,   1.48090574e-02,   9.77031394e-02,
                 8.03837049e-01,  4.55543734e+00,   2.13849673e+01]

synchrotron_comp = [6.79173120e-08,   2.32258906e-08,   6.64662670e-09,
                    1.54403227e-09,   4.19557011e-10,   1.07904119e-10]

CO_comp = [1*1.42*244.1, 0, 0.4*4.5*483.69,
           0.2*17.37*287.45, 0, 0]
