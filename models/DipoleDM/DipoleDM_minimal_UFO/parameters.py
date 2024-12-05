# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 12.1.0 for Linux x86 (64-bit) (March 18, 2020)
# Date: Thu 5 Dec 2024 11:40:29



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
LambdaUV = Parameter(name = 'LambdaUV',
                     nature = 'external',
                     type = 'real',
                     value = 5000.,
                     texname = '\\lambda',
                     lhablock = 'NPINPUTS',
                     lhacode = [ 1 ])

Caxx10 = Parameter(name = 'Caxx10',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{axx10}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 4 ])

Chxx10 = Parameter(name = 'Chxx10',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{hxx10}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 7 ])

ychi1 = Parameter(name = 'ychi1',
                  nature = 'external',
                  type = 'real',
                  value = 1.,
                  texname = 'y_{\\text{chi1}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 8 ])

sina = Parameter(name = 'sina',
                 nature = 'external',
                 type = 'real',
                 value = 0.2,
                 texname = '\\text{sina}',
                 lhablock = 'NPINPUTS',
                 lhacode = [ 11 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

vevD = Parameter(name = 'vevD',
                 nature = 'external',
                 type = 'real',
                 value = 1000.,
                 texname = '\\text{vevD}',
                 lhablock = 'FRBlock',
                 lhacode = [ 1 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

M0 = Parameter(name = 'M0',
               nature = 'external',
               type = 'real',
               value = 425.,
               texname = '\\text{M0}',
               lhablock = 'MASS',
               lhacode = [ 5000022 ])

M1 = Parameter(name = 'M1',
               nature = 'external',
               type = 'real',
               value = 500.,
               texname = '\\text{M1}',
               lhablock = 'MASS',
               lhacode = [ 5000023 ])

MSd = Parameter(name = 'MSd',
                nature = 'external',
                type = 'real',
                value = 1500,
                texname = '\\text{MSd}',
                lhablock = 'MASS',
                lhacode = [ 55 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.00407,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

w0 = Parameter(name = 'w0',
               nature = 'external',
               type = 'real',
               value = 0.,
               texname = '\\text{w0}',
               lhablock = 'DECAY',
               lhacode = [ 5000022 ])

w1 = Parameter(name = 'w1',
               nature = 'external',
               type = 'real',
               value = 0.00017,
               texname = '\\text{w1}',
               lhablock = 'DECAY',
               lhacode = [ 5000023 ])

WS = Parameter(name = 'WS',
               nature = 'external',
               type = 'real',
               value = 103.,
               texname = '\\text{WS}',
               lhablock = 'DECAY',
               lhacode = [ 55 ])

cosa = Parameter(name = 'cosa',
                 nature = 'internal',
                 type = 'real',
                 value = 'cmath.sqrt(1 - sina**2)',
                 texname = '\\text{cosa}')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

Caxx1x2 = Parameter(name = 'Caxx1x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx10',
                    texname = '\\text{Caxx1x2}')

Chxx1x2 = Parameter(name = 'Chxx1x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx10',
                    texname = '\\text{Chxx1x2}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

ychi2x2 = Parameter(name = 'ychi2x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi1',
                    texname = '\\text{ychi2x2}')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

lam2 = Parameter(name = 'lam2',
                 nature = 'internal',
                 type = 'real',
                 value = '(cosa**2*MSd**2)/(2.*vevD**2) + (MH**2*sina**2)/(2.*vevD**2)',
                 texname = '\\text{lam2}')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

GGH = Parameter(name = 'GGH',
                nature = 'internal',
                type = 'real',
                value = '-(G**2*cmath.sqrt(( ((6*MT**2)/MH**2 + (6*MT**2*cmath.asin(1/(2.*cmath.sqrt(MT**2/MH**2)))**2)/MH**2 - (24*MT**4*cmath.asin(1/(2.*cmath.sqrt(MT**2/MH**2)))**2)/MH**4)**2 if abs(MH**2/MT**2)/4.<1. else (36*MT**8*((cmath.pi**2*(-1 + MH**2/(4.*MT**2)) + MH**2/MT**2)**2 + 2*(cmath.pi**2*(-1 + MH**2/(4.*MT**2)) - MH**2/MT**2)*(-1 + MH**2/(4.*MT**2))*cmath.log((1 + 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2))/(1 - 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2)))**2 + (-1 + MH**2/(4.*MT**2))**2*cmath.log((1 + 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2))/(1 - 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2)))**4))/MH**8 )))/(12.*cmath.pi**2*vev)',
                texname = 'G_H')

lam1 = Parameter(name = 'lam1',
                 nature = 'internal',
                 type = 'real',
                 value = '(cosa**2*MH**2)/(2.*vev**2) + (MSd**2*sina**2)/(2.*vev**2)',
                 texname = '\\text{lam1}')

lam3 = Parameter(name = 'lam3',
                 nature = 'internal',
                 type = 'real',
                 value = '(cosa*(-MH**2 + MSd**2)*sina)/(vev*vevD)',
                 texname = '\\text{lam3}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/vev',
               texname = '\\text{yb}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

GGS = Parameter(name = 'GGS',
                nature = 'internal',
                type = 'real',
                value = '-(G**2*sina*cmath.sqrt(( ((6*MT**2)/MSd**2 + (6*MT**2*cmath.asin(1/(2.*cmath.sqrt(MT**2/MSd**2)))**2)/MSd**2 - (24*MT**4*cmath.asin(1/(2.*cmath.sqrt(MT**2/MSd**2)))**2)/MSd**4)**2 if abs(MSd**2/MT**2)/4.<1. else (36*MT**8*((cmath.pi**2*(-1 + MSd**2/(4.*MT**2)) + MSd**2/MT**2)**2 + 2*(cmath.pi**2*(-1 + MSd**2/(4.*MT**2)) - MSd**2/MT**2)*(-1 + MSd**2/(4.*MT**2))*cmath.log((1 + 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2))/(1 - 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2)))**2 + (-1 + MSd**2/(4.*MT**2))**2*cmath.log((1 + 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2))/(1 - 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2)))**4))/MSd**8 )))/(12.*cmath.pi**2*vev)',
                texname = 'G_S')

mu2h = Parameter(name = 'mu2h',
                 nature = 'internal',
                 type = 'real',
                 value = '-(lam1*vev**2) - (lam3*vevD**2)/2.',
                 texname = '\\mu _H')

mu2sd = Parameter(name = 'mu2sd',
                  nature = 'internal',
                  type = 'real',
                  value = '-(lam3*vev**2)/2. - lam2*vevD**2',
                  texname = '\\mu _2')

I1b33 = Parameter(name = 'I1b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I1b33}')

I2b33 = Parameter(name = 'I2b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I2b33}')

I3b33 = Parameter(name = 'I3b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I3b33}')

I4b33 = Parameter(name = 'I4b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I4b33}')

