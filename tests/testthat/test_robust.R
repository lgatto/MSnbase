context("Robust summarisation aggregation")

## This are the values calculated using data(msnset) and the orginal
## code provided by Adriaan Sticker, prior to refactoring to include
## into MSnbase. These data have been produced after removing X46 that
## contains one NA.
exp2 <- structure(c(885.188718025484, 17593.5475899093,
                    4923.62812560424, 1524.14767264389,
                    1069.94458515849, 1101.0620483011,
                    5852.6654376369, 760.032042354345,
                    6448.08286505193, 2278.3661184581,
                    32838.0439499617, 5135.5601051636,
                    21262.1475404352, 3715.08931949269,
                    4483.31990008056, 1147.0379543528,
                    26143.7541643903, 739.98608353734,
                    2770.37693640403, 2092.47630918864,
                    30047.2436312884, 11947.8807705827,
                    1259.1889827475, 4458.77171226405,
                    3824.35622986034, 2502.77917635441,
                    8635.31649124995, 3086.44980508555,
                    6670.69421123341, 31892.8927857503,
                    8504.48166594654, 2800.39559176937,
                    7540.63725898787, 1098.20952160703,
                    27638.3581933826, 1566.74773280695,
                    722.298177056015, 10154.9532804489,
                    11673.5920226627, 38262.1780306886,
                    1794.1744584286, 18545.6197157167,
                    5557.81802197453, 1399.89688059501,
                    1035.68883815967, 1124.16708243825,
                    5964.39464343898, 866.440367104486,
                    6234.19567901641, 2570.08000953496,
                    37066.0577703491, 5519.05329268798,
                    23168.7291266322, 4254.32256129384,
                    4873.99605610222, 1281.25856410339,
                    29677.4781315885, 799.350124068558,
                    2907.22208145075, 1993.87953759916,
                    31498.2176392823, 13061.8749812692,
                    1448.3734087348, 4786.51826088317,
                    3746.26739286631, 3013.04268418811,
                    10036.5290785991, 3266.41942225583,
                    7141.72540647816, 33634.69800926,
                    9924.07051063143, 3244.97753785551,
                    7679.15758221596, 984.645112269558,
                    33394.0252374262, 1444.27953283116,
                    1058.87115819752, 10486.9426673427,
                    11936.426138571, 32907.0889929717,
                    3157.01358505773, 19361.8367030621,
                    5775.20250520669, 1547.2184817791,
                    1029.41993487254, 1140.0931108892,
                    5970.36942510493, 883.909954763949,
                    6902.8902874887, 2785.56561793108,
                    41429.6267058477, 5828.37084272504,
                    25407.0684720129, 4748.46151560312,
                    6743.44063210487, 1278.70097099617,
                    29089.0593050867, 712.59829556942,
                    3055.51652965136, 2371.34243719932,
                    38489.8088214397, 12809.491407752,
                    1653.74895876646, 4350.78472106531,
                    4285.02115055174, 3017.53818337992,
                    9254.43168258667, 4093.80327095091,
                    7091.15089240577, 37674.727212131,
                    10058.2659687344, 2965.31877342984,
                    7561.98333513737, 1236.19285256579,
                    32104.2878604829, 1455.8965883255,
                    851.470159143209, 11018.1911818683,
                    12090.9319498833, 15000.1110553568,
                    7098.84112878367, 18328.2365038097,
                    5079.2951790709, 1563.22990202904,
                    999.695685420185, 1191.80548018217,
                    6006.8593209032, 791.329627946019,
                    6437.23029309139, 2446.79895107448,
                    39700.4751771986, 5583.72953182459,
                    25949.9540396184, 5249.90356134251,
                    4601.37757816911, 1175.08044352382,
                    27902.5607625172, 940.679303489625,
                    2946.79965814017, 2620.85662680864,
                    38275.2637470067, 12911.4785282575,
                    1261.52449379861, 4303.02562886477,
                    3649.44421585836, 2996.31272624433,
                    7769.74907445163, 3905.62405605242,
                    6779.99942643754, 37227.711932838,
                    9611.86445650645, 3184.19385924935,
                    8005.17654272541, 1173.12220480572,
                    26628.7277633697, 1557.19896062836,
                    944.628683254123, 11289.5520156138,
                    12268.6965850203, 10733.3352748144),
                  .Dim = c(40L, 4L),
                  .Dimnames = list(c("BSA", "ECA0172", "ECA0435",
                                     "ECA0452", "ECA0469", "ECA0621",
                                     "ECA0631", "ECA0691", "ECA0871",
                                     "ECA0978", "ECA1032", "ECA1093",
                                     "ECA1104", "ECA1294", "ECA1362",
                                     "ECA1363", "ECA1364", "ECA1422",
                                     "ECA1443", "ECA2186", "ECA2391",
                                     "ECA2421", "ECA2831", "ECA3082",
                                     "ECA3175", "ECA3349", "ECA3356",
                                     "ECA3377", "ECA3566", "ECA3882",
                                     "ECA3929", "ECA3969", "ECA4013",
                                     "ECA4026", "ECA4030", "ECA4037",
                                     "ECA4512", "ECA4513", "ECA4514",
                                     "ENO"),
                                   c("iTRAQ4.114", "iTRAQ4.115",
                                     "iTRAQ4.116", "iTRAQ4.117")))

## These have been created when only the single NA is removed
exp1 <- structure(c(885.188718025484, 17593.5475899093,
                    4923.62812560424, 1524.14767264389,
                    1069.94458515849, 1101.0620483011,
                    5852.6654376369, 760.032042354345,
                    6448.08286505193, 2278.3661184581,
                    32838.0439499617, 5135.5601051636,
                    21262.1475404352, 3715.08931949269,
                    4483.31990008056, 1147.0379543528,
                    26143.7541643903, 739.98608353734,
                    2770.37693640403, 2092.47630918864,
                    30047.2436312884, 11947.8807705827,
                    1259.1889827475, 4458.77171226405,
                    3824.35622986034, 2502.77917635441,
                    8635.31649124995, 3086.44980508555,
                    6670.69421123341, 31892.8927857503,
                    8504.48166594654, 2800.39559176937,
                    7540.63725898787, 1098.20952160703,
                    27638.3581933826, 1566.74773280695,
                    722.298177056015, 10154.9532804489,
                    11673.5920226627, 22671.7355697837,
                    1794.1744584286, 18545.6197157167,
                    5557.81802197453, 1399.89688059501,
                    1035.68883815967, 1124.16708243825,
                    5964.39464343898, 866.440367104486,
                    6234.19567901641, 2570.08000953496,
                    37066.0577703491, 5519.05329268798,
                    23168.7291266322, 4254.32256129384,
                    4873.99605610222, 1281.25856410339,
                    29677.4781315885, 799.350124068558,
                    2907.22208145075, 1993.87953759916,
                    31498.2176392823, 13061.8749812692,
                    1448.3734087348, 4786.51826088317,
                    3746.26739286631, 3013.04268418811,
                    10036.5290785991, 3266.41942225583,
                    7141.72540647816, 33634.69800926,
                    9924.07051063143, 3244.97753785551,
                    7679.15758221596, 984.645112269558,
                    33394.0252374262, 1444.27953283116,
                    1058.87115819752, 10486.9426673427,
                    11936.426138571, 19986.1636445553,
                    3157.01358505773, 19361.8367030621,
                    5775.20250520669, 1547.2184817791,
                    1029.41993487254, 1140.0931108892,
                    5970.36942510493, 883.909954763949,
                    6902.8902874887, 2785.56561793108,
                    41429.6267058477, 5828.37084272504,
                    25407.0684720129, 4748.46151560312,
                    6743.44063210487, 1278.70097099617,
                    29089.0593050867, 712.59829556942,
                    3055.51652965136, 2371.34243719932,
                    38489.8088214397, 12809.491407752,
                    1653.74895876646, 4350.78472106531,
                    4285.02115055174, 3017.53818337992,
                    9254.43168258667, 4093.80327095091,
                    7091.15089240577, 37674.727212131,
                    10058.2659687344, 2965.31877342984,
                    7561.98333513737, 1236.19285256579,
                    32104.2878604829, 1455.8965883255,
                    851.470159143209, 11018.1911818683,
                    12090.9319498833, 16679.5719334398,
                    7098.84112878367, 18328.2365038097,
                    5079.2951790709, 1563.22990202904,
                    999.695685420185, 1191.80548018217,
                    6006.8593209032, 791.329627946019,
                    6437.23029309139, 2446.79895107448,
                    39700.4751771986, 5583.72953182459,
                    25949.9540396184, 5249.90356134251,
                    4601.37757816911, 1175.08044352382,
                    27902.5607625172, 940.679303489625,
                    2946.79965814017, 2620.85662680864,
                    38275.2637470067, 12911.4785282575,
                    1261.52449379861, 4303.02562886477,
                    3649.44421585836, 2996.31272624433,
                    7769.74907445163, 3905.62405605242,
                    6779.99942643754, 37227.711932838,
                    9611.86445650645, 3184.19385924935,
                    8005.17654272541, 1173.12220480572,
                    26628.7277633697, 1557.19896062836,
                    944.628683254123, 11289.5520156138,
                    12268.6965850203, 16885.1437447593),
                  .Dim = c(40L, 4L),
                  .Dimnames = list(c("BSA", "ECA0172", "ECA0435",
                                     "ECA0452", "ECA0469", "ECA0621",
                                     "ECA0631", "ECA0691", "ECA0871",
                                     "ECA0978", "ECA1032", "ECA1093",
                                     "ECA1104", "ECA1294", "ECA1362",
                                     "ECA1363", "ECA1364", "ECA1422",
                                     "ECA1443", "ECA2186", "ECA2391",
                                     "ECA2421", "ECA2831", "ECA3082",
                                     "ECA3175", "ECA3349", "ECA3356",
                                     "ECA3377", "ECA3566", "ECA3882",
                                     "ECA3929", "ECA3969", "ECA4013",
                                     "ECA4026", "ECA4030", "ECA4037",
                                     "ECA4512", "ECA4513", "ECA4514",
                                     "ENO"),
                                   c("iTRAQ4.114", "iTRAQ4.115",
                                     "iTRAQ4.116", "iTRAQ4.117")))


test_that("Robust summarisation", {
    data("msnset")
    ## Single NA value is removed prior to fit
    expect_message(res1 <- combineFeatures(msnset,
                                           fcol = "ProteinAccession",
                                           method = "robust",
                                           maxit = 100L))
    expect_equal(exp1, exprs(res1))
    ## remove feature with missing value
    msnset <- filterNA(msnset)
    res2 <- combineFeatures(msnset,
                            fcol = "ProteinAccession",
                            method = "robust",
                            maxit = 100L)

    expect_equal(exp2, exprs(res2))

    notENO <- featureNames(res1) != "ENO"
    expect_equal(exprs(res1)[notENO, ], exprs(res2)[notENO, ])
    ## all false (0)
    expect_identical(sum(exprs(res1)[!notENO, ] == exprs(res2)[!notENO, ]), 0L)
})
