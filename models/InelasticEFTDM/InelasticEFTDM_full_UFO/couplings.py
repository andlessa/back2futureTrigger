# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 13.3.1 for Linux x86 (64-bit) (July 24, 2023)
# Date: Wed 23 Sep 2026 15:19:51


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-0.3333333333333333*(ee*complex(0,1))',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'ee**2*complex(0,1)',
                order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = '2*ee**2*complex(0,1)',
                order = {'QED':2})

GC_7 = Coupling(name = 'GC_7',
                value = '(ee**2*complex(0,1))/(2.*cw)',
                order = {'QED':2})

GC_8 = Coupling(name = 'GC_8',
                value = '-G',
                order = {'QCD':1})

GC_9 = Coupling(name = 'GC_9',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_10 = Coupling(name = 'GC_10',
                 value = 'complex(0,1)*G**2',
                 order = {'QCD':2})

GC_11 = Coupling(name = 'GC_11',
                 value = '-(complex(0,1)*GGH)',
                 order = {'QCD':2})

GC_12 = Coupling(name = 'GC_12',
                 value = '-(G*GGH)',
                 order = {'QCD':3})

GC_13 = Coupling(name = 'GC_13',
                 value = 'complex(0,1)*G**2*GGH',
                 order = {'QCD':4})

GC_14 = Coupling(name = 'GC_14',
                 value = '-(complex(0,1)*GGS)',
                 order = {'NP':1,'QCD':2})

GC_15 = Coupling(name = 'GC_15',
                 value = '-(G*GGS)',
                 order = {'NP':1,'QCD':3})

GC_16 = Coupling(name = 'GC_16',
                 value = 'complex(0,1)*G**2*GGS',
                 order = {'NP':1,'QCD':4})

GC_17 = Coupling(name = 'GC_17',
                 value = 'I1a33',
                 order = {'QED':1})

GC_18 = Coupling(name = 'GC_18',
                 value = '-I2a33',
                 order = {'QED':1})

GC_19 = Coupling(name = 'GC_19',
                 value = 'I3a33',
                 order = {'QED':1})

GC_20 = Coupling(name = 'GC_20',
                 value = '-I4a33',
                 order = {'QED':1})

GC_21 = Coupling(name = 'GC_21',
                 value = '-2*complex(0,1)*lam1',
                 order = {'QED':2})

GC_22 = Coupling(name = 'GC_22',
                 value = '-4*complex(0,1)*lam1',
                 order = {'QED':2})

GC_23 = Coupling(name = 'GC_23',
                 value = '-6*complex(0,1)*lam1',
                 order = {'QED':2})

GC_24 = Coupling(name = 'GC_24',
                 value = '-6*complex(0,1)*lam2',
                 order = {'QED':2})

GC_25 = Coupling(name = 'GC_25',
                 value = '-(complex(0,1)*lam3)',
                 order = {'QED':2})

GC_26 = Coupling(name = 'GC_26',
                 value = '-((Chxx00*complex(0,1))/LambdaUV)',
                 order = {'NP':1})

GC_27 = Coupling(name = 'GC_27',
                 value = '-((Chxx10*complex(0,1))/LambdaUV)',
                 order = {'NP':1})

GC_28 = Coupling(name = 'GC_28',
                 value = '-((Chxx11*complex(0,1))/LambdaUV)',
                 order = {'NP':1})

GC_29 = Coupling(name = 'GC_29',
                 value = '-((Chxx20*complex(0,1))/LambdaUV)',
                 order = {'NP':1})

GC_30 = Coupling(name = 'GC_30',
                 value = '-((Chxx21*complex(0,1))/LambdaUV)',
                 order = {'NP':1})

GC_31 = Coupling(name = 'GC_31',
                 value = '-((Chxx22*complex(0,1))/LambdaUV)',
                 order = {'NP':1})

GC_32 = Coupling(name = 'GC_32',
                 value = '-0.5*(ee**2*sina)/cw',
                 order = {'NP':1,'QED':2})

GC_33 = Coupling(name = 'GC_33',
                 value = '(ee**2*sina)/(2.*cw)',
                 order = {'NP':1,'QED':2})

GC_34 = Coupling(name = 'GC_34',
                 value = '-((Chxx00*complex(0,1)*sina**2)/LambdaUV)',
                 order = {'NP':3})

GC_35 = Coupling(name = 'GC_35',
                 value = '(Chxx00*complex(0,1)*sina**2)/LambdaUV',
                 order = {'NP':3})

GC_36 = Coupling(name = 'GC_36',
                 value = '-((Chxx10*complex(0,1)*sina**2)/LambdaUV)',
                 order = {'NP':3})

GC_37 = Coupling(name = 'GC_37',
                 value = '(Chxx10*complex(0,1)*sina**2)/LambdaUV',
                 order = {'NP':3})

GC_38 = Coupling(name = 'GC_38',
                 value = '-((Chxx11*complex(0,1)*sina**2)/LambdaUV)',
                 order = {'NP':3})

GC_39 = Coupling(name = 'GC_39',
                 value = '(Chxx11*complex(0,1)*sina**2)/LambdaUV',
                 order = {'NP':3})

GC_40 = Coupling(name = 'GC_40',
                 value = '-((Chxx20*complex(0,1)*sina**2)/LambdaUV)',
                 order = {'NP':3})

GC_41 = Coupling(name = 'GC_41',
                 value = '(Chxx20*complex(0,1)*sina**2)/LambdaUV',
                 order = {'NP':3})

GC_42 = Coupling(name = 'GC_42',
                 value = '-((Chxx21*complex(0,1)*sina**2)/LambdaUV)',
                 order = {'NP':3})

GC_43 = Coupling(name = 'GC_43',
                 value = '(Chxx21*complex(0,1)*sina**2)/LambdaUV',
                 order = {'NP':3})

GC_44 = Coupling(name = 'GC_44',
                 value = '-((Chxx22*complex(0,1)*sina**2)/LambdaUV)',
                 order = {'NP':3})

GC_45 = Coupling(name = 'GC_45',
                 value = '(Chxx22*complex(0,1)*sina**2)/LambdaUV',
                 order = {'NP':3})

GC_46 = Coupling(name = 'GC_46',
                 value = '-0.5*(ee**2*cmath.sqrt(1 - sina**2))/cw',
                 order = {'QED':2})

GC_47 = Coupling(name = 'GC_47',
                 value = '(ee**2*cmath.sqrt(1 - sina**2))/(2.*cw)',
                 order = {'QED':2})

GC_48 = Coupling(name = 'GC_48',
                 value = '-((Chxx00*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/LambdaUV)',
                 order = {'NP':2})

GC_49 = Coupling(name = 'GC_49',
                 value = '-((Chxx10*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/LambdaUV)',
                 order = {'NP':2})

GC_50 = Coupling(name = 'GC_50',
                 value = '-((Chxx11*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/LambdaUV)',
                 order = {'NP':2})

GC_51 = Coupling(name = 'GC_51',
                 value = '-((Chxx20*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/LambdaUV)',
                 order = {'NP':2})

GC_52 = Coupling(name = 'GC_52',
                 value = '-((Chxx21*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/LambdaUV)',
                 order = {'NP':2})

GC_53 = Coupling(name = 'GC_53',
                 value = '-((Chxx22*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/LambdaUV)',
                 order = {'NP':2})

GC_54 = Coupling(name = 'GC_54',
                 value = '2*complex(0,1)*lam1*sina**2 - complex(0,1)*lam3*sina**2',
                 order = {'NP':2,'QED':2})

GC_55 = Coupling(name = 'GC_55',
                 value = '-2*complex(0,1)*lam1*sina**2 + complex(0,1)*lam3*sina**2',
                 order = {'NP':2,'QED':2})

GC_56 = Coupling(name = 'GC_56',
                 value = '12*complex(0,1)*lam1*sina**2 - 6*complex(0,1)*lam3*sina**2',
                 order = {'NP':2,'QED':2})

GC_57 = Coupling(name = 'GC_57',
                 value = '12*complex(0,1)*lam2*sina**2 - 6*complex(0,1)*lam3*sina**2',
                 order = {'NP':2,'QED':2})

GC_58 = Coupling(name = 'GC_58',
                 value = '-6*complex(0,1)*lam1*sina**2 - 6*complex(0,1)*lam2*sina**2 + 6*complex(0,1)*lam3*sina**2',
                 order = {'NP':2,'QED':2})

GC_59 = Coupling(name = 'GC_59',
                 value = '6*complex(0,1)*lam1*sina**4 + 6*complex(0,1)*lam2*sina**4 - 6*complex(0,1)*lam3*sina**4',
                 order = {'NP':4,'QED':2})

GC_60 = Coupling(name = 'GC_60',
                 value = '-6*complex(0,1)*lam1*sina**4 - 6*complex(0,1)*lam2*sina**4 + 6*complex(0,1)*lam3*sina**4',
                 order = {'NP':4,'QED':2})

GC_61 = Coupling(name = 'GC_61',
                 value = '-2*complex(0,1)*lam1*sina*cmath.sqrt(1 - sina**2) + complex(0,1)*lam3*sina*cmath.sqrt(1 - sina**2)',
                 order = {'NP':1,'QED':2})

GC_62 = Coupling(name = 'GC_62',
                 value = '6*complex(0,1)*lam2*sina*cmath.sqrt(1 - sina**2) - 3*complex(0,1)*lam3*sina*cmath.sqrt(1 - sina**2)',
                 order = {'NP':1,'QED':2})

GC_63 = Coupling(name = 'GC_63',
                 value = '-6*complex(0,1)*lam1*sina*cmath.sqrt(1 - sina**2) + 3*complex(0,1)*lam3*sina*cmath.sqrt(1 - sina**2)',
                 order = {'NP':1,'QED':2})

GC_64 = Coupling(name = 'GC_64',
                 value = '6*complex(0,1)*lam1*sina**3*cmath.sqrt(1 - sina**2) + 6*complex(0,1)*lam2*sina**3*cmath.sqrt(1 - sina**2) - 6*complex(0,1)*lam3*sina**3*cmath.sqrt(1 - sina**2)',
                 order = {'NP':3,'QED':2})

GC_65 = Coupling(name = 'GC_65',
                 value = '-6*complex(0,1)*lam1*sina**3*cmath.sqrt(1 - sina**2) - 6*complex(0,1)*lam2*sina**3*cmath.sqrt(1 - sina**2) + 6*complex(0,1)*lam3*sina**3*cmath.sqrt(1 - sina**2)',
                 order = {'NP':3,'QED':2})

GC_66 = Coupling(name = 'GC_66',
                 value = '(ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_67 = Coupling(name = 'GC_67',
                 value = '-((ee**2*complex(0,1))/sw**2)',
                 order = {'QED':2})

GC_68 = Coupling(name = 'GC_68',
                 value = '(cw**2*ee**2*complex(0,1))/sw**2',
                 order = {'QED':2})

GC_69 = Coupling(name = 'GC_69',
                 value = '(-2*cw**2*ee**2*complex(0,1))/sw**2',
                 order = {'QED':2})

GC_70 = Coupling(name = 'GC_70',
                 value = '-0.5*(ee**2*complex(0,1)*sina**2)/sw**2',
                 order = {'NP':2,'QED':2})

GC_71 = Coupling(name = 'GC_71',
                 value = '(ee**2*complex(0,1)*sina**2)/(2.*sw**2)',
                 order = {'NP':2,'QED':2})

GC_72 = Coupling(name = 'GC_72',
                 value = '(ee**2*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/(2.*sw**2)',
                 order = {'NP':1,'QED':2})

GC_73 = Coupling(name = 'GC_73',
                 value = '-0.5*(ee*complex(0,1))/sw',
                 order = {'QED':1})

GC_74 = Coupling(name = 'GC_74',
                 value = '(ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_75 = Coupling(name = 'GC_75',
                 value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_76 = Coupling(name = 'GC_76',
                 value = '-((cw*ee*complex(0,1))/sw)',
                 order = {'QED':1})

GC_77 = Coupling(name = 'GC_77',
                 value = '(cw*ee*complex(0,1))/sw',
                 order = {'QED':1})

GC_78 = Coupling(name = 'GC_78',
                 value = '-0.5*(ee**2*complex(0,1))/sw',
                 order = {'QED':2})

GC_79 = Coupling(name = 'GC_79',
                 value = '(cw*ee**2*complex(0,1))/sw',
                 order = {'QED':2})

GC_80 = Coupling(name = 'GC_80',
                 value = '(-2*cw*ee**2*complex(0,1))/sw',
                 order = {'QED':2})

GC_81 = Coupling(name = 'GC_81',
                 value = '-0.5*(ee*sina)/sw',
                 order = {'NP':1,'QED':1})

GC_82 = Coupling(name = 'GC_82',
                 value = '-0.5*(ee**2*sina)/sw',
                 order = {'NP':1,'QED':2})

GC_83 = Coupling(name = 'GC_83',
                 value = '(ee**2*sina)/(2.*sw)',
                 order = {'NP':1,'QED':2})

GC_84 = Coupling(name = 'GC_84',
                 value = '-0.5*(ee*cmath.sqrt(1 - sina**2))/sw',
                 order = {'QED':1})

GC_85 = Coupling(name = 'GC_85',
                 value = '-0.5*(ee**2*cmath.sqrt(1 - sina**2))/sw',
                 order = {'QED':2})

GC_86 = Coupling(name = 'GC_86',
                 value = '(ee**2*cmath.sqrt(1 - sina**2))/(2.*sw)',
                 order = {'QED':2})

GC_87 = Coupling(name = 'GC_87',
                 value = '(ee*complex(0,1)*sw)/(3.*cw)',
                 order = {'QED':1})

GC_88 = Coupling(name = 'GC_88',
                 value = '(-2*ee*complex(0,1)*sw)/(3.*cw)',
                 order = {'QED':1})

GC_89 = Coupling(name = 'GC_89',
                 value = '(ee*complex(0,1)*sw)/cw',
                 order = {'QED':1})

GC_90 = Coupling(name = 'GC_90',
                 value = '-0.5*(cw*ee*complex(0,1))/sw - (ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_91 = Coupling(name = 'GC_91',
                 value = '(cw*ee*complex(0,1))/(2.*sw) - (ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_92 = Coupling(name = 'GC_92',
                 value = '-0.5*(cw*ee*complex(0,1))/sw + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_93 = Coupling(name = 'GC_93',
                 value = '(cw*ee*complex(0,1))/(2.*sw) + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_94 = Coupling(name = 'GC_94',
                 value = '(cw*ee**2*complex(0,1))/sw - (ee**2*complex(0,1)*sw)/cw',
                 order = {'QED':2})

GC_95 = Coupling(name = 'GC_95',
                 value = '-0.5*(cw*ee*sina)/sw - (ee*sina*sw)/(2.*cw)',
                 order = {'NP':1,'QED':1})

GC_96 = Coupling(name = 'GC_96',
                 value = '-0.5*(cw*ee*cmath.sqrt(1 - sina**2))/sw - (ee*sw*cmath.sqrt(1 - sina**2))/(2.*cw)',
                 order = {'QED':1})

GC_97 = Coupling(name = 'GC_97',
                 value = '-(ee**2*complex(0,1)) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_98 = Coupling(name = 'GC_98',
                 value = 'ee**2*complex(0,1) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_99 = Coupling(name = 'GC_99',
                 value = '-(ee**2*complex(0,1)*sina**2) - (cw**2*ee**2*complex(0,1)*sina**2)/(2.*sw**2) - (ee**2*complex(0,1)*sina**2*sw**2)/(2.*cw**2)',
                 order = {'NP':2,'QED':2})

GC_100 = Coupling(name = 'GC_100',
                  value = 'ee**2*complex(0,1)*sina**2 + (cw**2*ee**2*complex(0,1)*sina**2)/(2.*sw**2) + (ee**2*complex(0,1)*sina**2*sw**2)/(2.*cw**2)',
                  order = {'NP':2,'QED':2})

GC_101 = Coupling(name = 'GC_101',
                  value = 'ee**2*complex(0,1)*sina*cmath.sqrt(1 - sina**2) + (cw**2*ee**2*complex(0,1)*sina*cmath.sqrt(1 - sina**2))/(2.*sw**2) + (ee**2*complex(0,1)*sina*sw**2*cmath.sqrt(1 - sina**2))/(2.*cw**2)',
                  order = {'NP':1,'QED':2})

GC_102 = Coupling(name = 'GC_102',
                  value = '-0.5*(ee**2*vev)/cw',
                  order = {'QED':1})

GC_103 = Coupling(name = 'GC_103',
                  value = '(ee**2*vev)/(2.*cw)',
                  order = {'QED':1})

GC_104 = Coupling(name = 'GC_104',
                  value = '-2*complex(0,1)*lam1*sina*vev',
                  order = {'NP':1,'QED':1})

GC_105 = Coupling(name = 'GC_105',
                  value = '-3*complex(0,1)*lam3*sina*vev',
                  order = {'NP':1,'QED':1})

GC_106 = Coupling(name = 'GC_106',
                  value = '-((Chxx00*complex(0,1)*sina*vev)/LambdaUV)',
                  order = {'NP':2,'QED':-1})

GC_107 = Coupling(name = 'GC_107',
                  value = '-((Chxx10*complex(0,1)*sina*vev)/LambdaUV)',
                  order = {'NP':2,'QED':-1})

GC_108 = Coupling(name = 'GC_108',
                  value = '-((Chxx11*complex(0,1)*sina*vev)/LambdaUV)',
                  order = {'NP':2,'QED':-1})

GC_109 = Coupling(name = 'GC_109',
                  value = '-((Chxx20*complex(0,1)*sina*vev)/LambdaUV)',
                  order = {'NP':2,'QED':-1})

GC_110 = Coupling(name = 'GC_110',
                  value = '-((Chxx21*complex(0,1)*sina*vev)/LambdaUV)',
                  order = {'NP':2,'QED':-1})

GC_111 = Coupling(name = 'GC_111',
                  value = '-((Chxx22*complex(0,1)*sina*vev)/LambdaUV)',
                  order = {'NP':2,'QED':-1})

GC_112 = Coupling(name = 'GC_112',
                  value = '-2*complex(0,1)*lam1*vev*cmath.sqrt(1 - sina**2)',
                  order = {'QED':1})

GC_113 = Coupling(name = 'GC_113',
                  value = '-6*complex(0,1)*lam1*vev*cmath.sqrt(1 - sina**2)',
                  order = {'QED':1})

GC_114 = Coupling(name = 'GC_114',
                  value = '-(complex(0,1)*lam3*vev*cmath.sqrt(1 - sina**2))',
                  order = {'QED':1})

GC_115 = Coupling(name = 'GC_115',
                  value = '-((Chxx00*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/LambdaUV)',
                  order = {'NP':1,'QED':-1})

GC_116 = Coupling(name = 'GC_116',
                  value = '-((Chxx10*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/LambdaUV)',
                  order = {'NP':1,'QED':-1})

GC_117 = Coupling(name = 'GC_117',
                  value = '-((Chxx11*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/LambdaUV)',
                  order = {'NP':1,'QED':-1})

GC_118 = Coupling(name = 'GC_118',
                  value = '-((Chxx20*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/LambdaUV)',
                  order = {'NP':1,'QED':-1})

GC_119 = Coupling(name = 'GC_119',
                  value = '-((Chxx21*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/LambdaUV)',
                  order = {'NP':1,'QED':-1})

GC_120 = Coupling(name = 'GC_120',
                  value = '-((Chxx22*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/LambdaUV)',
                  order = {'NP':1,'QED':-1})

GC_121 = Coupling(name = 'GC_121',
                  value = '-0.25*(ee**2*vev)/sw**2',
                  order = {'QED':1})

GC_122 = Coupling(name = 'GC_122',
                  value = '(ee**2*vev)/(4.*sw**2)',
                  order = {'QED':1})

GC_123 = Coupling(name = 'GC_123',
                  value = '-0.25*(ee**2*complex(0,1)*sina*vev)/sw**2',
                  order = {'NP':1,'QED':1})

GC_124 = Coupling(name = 'GC_124',
                  value = '(ee**2*complex(0,1)*sina*vev)/(2.*sw**2)',
                  order = {'NP':1,'QED':1})

GC_125 = Coupling(name = 'GC_125',
                  value = '-0.25*(ee**2*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/sw**2',
                  order = {'QED':1})

GC_126 = Coupling(name = 'GC_126',
                  value = '(ee**2*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/(2.*sw**2)',
                  order = {'QED':1})

GC_127 = Coupling(name = 'GC_127',
                  value = '-0.5*(ee**2*vev)/sw',
                  order = {'QED':1})

GC_128 = Coupling(name = 'GC_128',
                  value = '(ee**2*vev)/(2.*sw)',
                  order = {'QED':1})

GC_129 = Coupling(name = 'GC_129',
                  value = '-6*complex(0,1)*lam1*sina*vev + 2*complex(0,1)*lam3*sina*vev',
                  order = {'NP':1,'QED':1})

GC_130 = Coupling(name = 'GC_130',
                  value = '6*complex(0,1)*lam1*sina**3*vev - 3*complex(0,1)*lam3*sina**3*vev',
                  order = {'NP':3,'QED':1})

GC_131 = Coupling(name = 'GC_131',
                  value = '-6*complex(0,1)*lam1*sina**3*vev + 3*complex(0,1)*lam3*sina**3*vev',
                  order = {'NP':3,'QED':1})

GC_132 = Coupling(name = 'GC_132',
                  value = '6*complex(0,1)*lam1*sina**2*vev*cmath.sqrt(1 - sina**2) - 3*complex(0,1)*lam3*sina**2*vev*cmath.sqrt(1 - sina**2)',
                  order = {'NP':2,'QED':1})

GC_133 = Coupling(name = 'GC_133',
                  value = '-6*complex(0,1)*lam1*sina**2*vev*cmath.sqrt(1 - sina**2) + 3*complex(0,1)*lam3*sina**2*vev*cmath.sqrt(1 - sina**2)',
                  order = {'NP':2,'QED':1})

GC_134 = Coupling(name = 'GC_134',
                  value = '-0.25*(ee**2*vev)/cw - (cw*ee**2*vev)/(4.*sw**2)',
                  order = {'QED':1})

GC_135 = Coupling(name = 'GC_135',
                  value = '(ee**2*vev)/(4.*cw) - (cw*ee**2*vev)/(4.*sw**2)',
                  order = {'QED':1})

GC_136 = Coupling(name = 'GC_136',
                  value = '-0.25*(ee**2*vev)/cw + (cw*ee**2*vev)/(4.*sw**2)',
                  order = {'QED':1})

GC_137 = Coupling(name = 'GC_137',
                  value = '(ee**2*vev)/(4.*cw) + (cw*ee**2*vev)/(4.*sw**2)',
                  order = {'QED':1})

GC_138 = Coupling(name = 'GC_138',
                  value = '-0.5*(ee**2*complex(0,1)*sina*vev) - (cw**2*ee**2*complex(0,1)*sina*vev)/(4.*sw**2) - (ee**2*complex(0,1)*sina*sw**2*vev)/(4.*cw**2)',
                  order = {'NP':1,'QED':1})

GC_139 = Coupling(name = 'GC_139',
                  value = 'ee**2*complex(0,1)*sina*vev + (cw**2*ee**2*complex(0,1)*sina*vev)/(2.*sw**2) + (ee**2*complex(0,1)*sina*sw**2*vev)/(2.*cw**2)',
                  order = {'NP':1,'QED':1})

GC_140 = Coupling(name = 'GC_140',
                  value = '-0.5*(ee**2*complex(0,1)*vev*cmath.sqrt(1 - sina**2)) - (cw**2*ee**2*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/(4.*sw**2) - (ee**2*complex(0,1)*sw**2*vev*cmath.sqrt(1 - sina**2))/(4.*cw**2)',
                  order = {'QED':1})

GC_141 = Coupling(name = 'GC_141',
                  value = 'ee**2*complex(0,1)*vev*cmath.sqrt(1 - sina**2) + (cw**2*ee**2*complex(0,1)*vev*cmath.sqrt(1 - sina**2))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2*vev*cmath.sqrt(1 - sina**2))/(2.*cw**2)',
                  order = {'QED':1})

GC_142 = Coupling(name = 'GC_142',
                  value = 'complex(0,1)*lam3*sina*vevD',
                  order = {'NP':1,'QED':2})

GC_143 = Coupling(name = 'GC_143',
                  value = '3*complex(0,1)*lam3*sina*vevD',
                  order = {'NP':1,'QED':2})

GC_144 = Coupling(name = 'GC_144',
                  value = '-6*complex(0,1)*lam2*vevD*cmath.sqrt(1 - sina**2)',
                  order = {'QED':2})

GC_145 = Coupling(name = 'GC_145',
                  value = '-(complex(0,1)*lam3*vevD*cmath.sqrt(1 - sina**2))',
                  order = {'QED':2})

GC_146 = Coupling(name = 'GC_146',
                  value = '6*complex(0,1)*lam2*sina*vevD - 2*complex(0,1)*lam3*sina*vevD',
                  order = {'NP':1,'QED':2})

GC_147 = Coupling(name = 'GC_147',
                  value = '6*complex(0,1)*lam2*sina**3*vevD - 3*complex(0,1)*lam3*sina**3*vevD',
                  order = {'NP':3,'QED':2})

GC_148 = Coupling(name = 'GC_148',
                  value = '-6*complex(0,1)*lam2*sina**3*vevD + 3*complex(0,1)*lam3*sina**3*vevD',
                  order = {'NP':3,'QED':2})

GC_149 = Coupling(name = 'GC_149',
                  value = '6*complex(0,1)*lam2*sina**2*vevD*cmath.sqrt(1 - sina**2) - 3*complex(0,1)*lam3*sina**2*vevD*cmath.sqrt(1 - sina**2)',
                  order = {'NP':2,'QED':2})

GC_150 = Coupling(name = 'GC_150',
                  value = '-6*complex(0,1)*lam2*sina**2*vevD*cmath.sqrt(1 - sina**2) + 3*complex(0,1)*lam3*sina**2*vevD*cmath.sqrt(1 - sina**2)',
                  order = {'NP':2,'QED':2})

GC_151 = Coupling(name = 'GC_151',
                  value = '-(yb/cmath.sqrt(2))',
                  order = {'QED':1})

GC_152 = Coupling(name = 'GC_152',
                  value = '-((complex(0,1)*sina*yb)/cmath.sqrt(2))',
                  order = {'NP':1,'QED':1})

GC_153 = Coupling(name = 'GC_153',
                  value = '-((complex(0,1)*yb*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'QED':1})

GC_154 = Coupling(name = 'GC_154',
                  value = '(complex(0,1)*sina*ychi00)/cmath.sqrt(2)',
                  order = {'NP':2})

GC_155 = Coupling(name = 'GC_155',
                  value = '-((complex(0,1)*ychi00*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'NP':1})

GC_156 = Coupling(name = 'GC_156',
                  value = '(complex(0,1)*sina*ychi10)/cmath.sqrt(2)',
                  order = {'NP':2})

GC_157 = Coupling(name = 'GC_157',
                  value = '-((complex(0,1)*ychi10*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'NP':1})

GC_158 = Coupling(name = 'GC_158',
                  value = '(complex(0,1)*sina*ychi11)/cmath.sqrt(2)',
                  order = {'NP':2})

GC_159 = Coupling(name = 'GC_159',
                  value = '-((complex(0,1)*ychi11*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'NP':1})

GC_160 = Coupling(name = 'GC_160',
                  value = '(complex(0,1)*sina*ychi20)/cmath.sqrt(2)',
                  order = {'NP':2})

GC_161 = Coupling(name = 'GC_161',
                  value = '-((complex(0,1)*ychi20*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'NP':1})

GC_162 = Coupling(name = 'GC_162',
                  value = '(complex(0,1)*sina*ychi21)/cmath.sqrt(2)',
                  order = {'NP':2})

GC_163 = Coupling(name = 'GC_163',
                  value = '-((complex(0,1)*ychi21*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'NP':1})

GC_164 = Coupling(name = 'GC_164',
                  value = '(complex(0,1)*sina*ychi22)/cmath.sqrt(2)',
                  order = {'NP':2})

GC_165 = Coupling(name = 'GC_165',
                  value = '-((complex(0,1)*ychi22*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'NP':1})

GC_166 = Coupling(name = 'GC_166',
                  value = 'yt/cmath.sqrt(2)',
                  order = {'QED':1})

GC_167 = Coupling(name = 'GC_167',
                  value = '-((complex(0,1)*sina*yt)/cmath.sqrt(2))',
                  order = {'NP':1,'QED':1})

GC_168 = Coupling(name = 'GC_168',
                  value = '-((complex(0,1)*yt*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'QED':1})

GC_169 = Coupling(name = 'GC_169',
                  value = '-ytau',
                  order = {'QED':1})

GC_170 = Coupling(name = 'GC_170',
                  value = 'ytau',
                  order = {'QED':1})

GC_171 = Coupling(name = 'GC_171',
                  value = '-(ytau/cmath.sqrt(2))',
                  order = {'QED':1})

GC_172 = Coupling(name = 'GC_172',
                  value = '-((complex(0,1)*sina*ytau)/cmath.sqrt(2))',
                  order = {'NP':1,'QED':1})

GC_173 = Coupling(name = 'GC_173',
                  value = '-((complex(0,1)*ytau*cmath.sqrt(1 - sina**2))/cmath.sqrt(2))',
                  order = {'QED':1})

