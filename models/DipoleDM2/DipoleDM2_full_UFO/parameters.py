# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 14.2.0 for Linux x86 (64-bit) (December 26, 2024)
# Date: Wed 13 Aug 2025 12:04:20



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
                     texname = '\\Lambda',
                     lhablock = 'NPINPUTS',
                     lhacode = [ 1 ])

Caxx2 = Parameter(name = 'Caxx2',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'C_{\\text{axx2}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 2 ])

Caxx1 = Parameter(name = 'Caxx1',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'C_{\\text{axx1}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 3 ])

Caxx0 = Parameter(name = 'Caxx0',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'C_{\\text{axx0}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 4 ])

Caxx10 = Parameter(name = 'Caxx10',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{axx10}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 5 ])

Caxx20 = Parameter(name = 'Caxx20',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{axx20}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 6 ])

Caxx21 = Parameter(name = 'Caxx21',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{axx21}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 7 ])

Chxx2 = Parameter(name = 'Chxx2',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'C_{\\text{hxx2}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 8 ])

Chxx1 = Parameter(name = 'Chxx1',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'C_{\\text{hxx1}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 9 ])

Chxx0 = Parameter(name = 'Chxx0',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'C_{\\text{hxx0}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 10 ])

Chxx10 = Parameter(name = 'Chxx10',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{hxx10}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 11 ])

Chxx20 = Parameter(name = 'Chxx20',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{hxx20}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 12 ])

Chxx21 = Parameter(name = 'Chxx21',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = 'C_{\\text{hxx21}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 13 ])

ychi2 = Parameter(name = 'ychi2',
                  nature = 'external',
                  type = 'real',
                  value = 1.,
                  texname = 'y_{\\text{chi2}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 14 ])

ychi1 = Parameter(name = 'ychi1',
                  nature = 'external',
                  type = 'real',
                  value = 1.,
                  texname = 'y_{\\text{chi1}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 15 ])

ychi0 = Parameter(name = 'ychi0',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'y_{\\text{chi0}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 16 ])

ychi10 = Parameter(name = 'ychi10',
                   nature = 'external',
                   type = 'real',
                   value = 0.,
                   texname = 'y_{\\text{chi10}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 17 ])

ychi20 = Parameter(name = 'ychi20',
                   nature = 'external',
                   type = 'real',
                   value = 0.,
                   texname = 'y_{\\text{chi20}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 18 ])

ychi21 = Parameter(name = 'ychi21',
                   nature = 'external',
                   type = 'real',
                   value = 0.,
                   texname = 'y_{\\text{chi21}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 19 ])

sina = Parameter(name = 'sina',
                 nature = 'external',
                 type = 'real',
                 value = 0.2,
                 texname = '\\text{sina}',
                 lhablock = 'NPINPUTS',
                 lhacode = [ 20 ])

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
               lhacode = [ 4000022 ])

M1 = Parameter(name = 'M1',
               nature = 'external',
               type = 'real',
               value = 500.,
               texname = '\\text{M1}',
               lhablock = 'MASS',
               lhacode = [ 4000023 ])

M2 = Parameter(name = 'M2',
               nature = 'external',
               type = 'real',
               value = 510.,
               texname = '\\text{M2}',
               lhablock = 'MASS',
               lhacode = [ 4000024 ])

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
               lhacode = [ 4000022 ])

w1 = Parameter(name = 'w1',
               nature = 'external',
               type = 'real',
               value = 0.00017,
               texname = '\\text{w1}',
               lhablock = 'DECAY',
               lhacode = [ 4000023 ])

w2 = Parameter(name = 'w2',
               nature = 'external',
               type = 'real',
               value = 0.00017,
               texname = '\\text{w2}',
               lhablock = 'DECAY',
               lhacode = [ 4000024 ])

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

Caxx1x1 = Parameter(name = 'Caxx1x1',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx0',
                    texname = '\\text{Caxx1x1}')

Caxx1x2 = Parameter(name = 'Caxx1x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx10',
                    texname = '\\text{Caxx1x2}')

Caxx1x3 = Parameter(name = 'Caxx1x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx20',
                    texname = '\\text{Caxx1x3}')

Caxx2x2 = Parameter(name = 'Caxx2x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx1',
                    texname = '\\text{Caxx2x2}')

Caxx2x3 = Parameter(name = 'Caxx2x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx21',
                    texname = '\\text{Caxx2x3}')

Caxx3x3 = Parameter(name = 'Caxx3x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'Caxx2',
                    texname = '\\text{Caxx3x3}')

Chxx1x1 = Parameter(name = 'Chxx1x1',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx0',
                    texname = '\\text{Chxx1x1}')

Chxx1x2 = Parameter(name = 'Chxx1x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx10',
                    texname = '\\text{Chxx1x2}')

Chxx1x3 = Parameter(name = 'Chxx1x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx20',
                    texname = '\\text{Chxx1x3}')

Chxx2x2 = Parameter(name = 'Chxx2x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx1',
                    texname = '\\text{Chxx2x2}')

Chxx2x3 = Parameter(name = 'Chxx2x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx21',
                    texname = '\\text{Chxx2x3}')

Chxx3x3 = Parameter(name = 'Chxx3x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'Chxx2',
                    texname = '\\text{Chxx3x3}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

ychi1x1 = Parameter(name = 'ychi1x1',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi0',
                    texname = '\\text{ychi1x1}')

ychi1x2 = Parameter(name = 'ychi1x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi10',
                    texname = '\\text{ychi1x2}')

ychi1x3 = Parameter(name = 'ychi1x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi20',
                    texname = '\\text{ychi1x3}')

ychi2x2 = Parameter(name = 'ychi2x2',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi1',
                    texname = '\\text{ychi2x2}')

ychi2x3 = Parameter(name = 'ychi2x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi21',
                    texname = '\\text{ychi2x3}')

ychi3x3 = Parameter(name = 'ychi3x3',
                    nature = 'internal',
                    type = 'real',
                    value = 'ychi2',
                    texname = '\\text{ychi3x3}')

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
                value = '-0.08333333333333333*(cosa*G**2*cmath.sqrt(( ((6*MT**2)/MH**2 + (6*MT**2*cmath.asin(1/(2.*cmath.sqrt(MT**2/MH**2)))**2)/MH**2 - (24*MT**4*cmath.asin(1/(2.*cmath.sqrt(MT**2/MH**2)))**2)/MH**4)**2 if abs(MH**2/MT**2)/4.<1. else (36*MT**8*((cmath.pi**2*(-1 + MH**2/(4.*MT**2)) + MH**2/MT**2)**2 + 2*(cmath.pi**2*(-1 + MH**2/(4.*MT**2)) - MH**2/MT**2)*(-1 + MH**2/(4.*MT**2))*cmath.log((1 + 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2))/(1 - 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2)))**2 + (-1 + MH**2/(4.*MT**2))**2*cmath.log((1 + 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2))/(1 - 2*cmath.sqrt(((-1 + MH**2/(4.*MT**2))*MT**2)/MH**2)))**4))/MH**8 )))/(cmath.pi**2*vev)',
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
                value = '-0.08333333333333333*(G**2*sina*cmath.sqrt(( ((6*MT**2)/MSd**2 + (6*MT**2*cmath.asin(1/(2.*cmath.sqrt(MT**2/MSd**2)))**2)/MSd**2 - (24*MT**4*cmath.asin(1/(2.*cmath.sqrt(MT**2/MSd**2)))**2)/MSd**4)**2 if abs(MSd**2/MT**2)/4.<1. else (36*MT**8*((cmath.pi**2*(-1 + MSd**2/(4.*MT**2)) + MSd**2/MT**2)**2 + 2*(cmath.pi**2*(-1 + MSd**2/(4.*MT**2)) - MSd**2/MT**2)*(-1 + MSd**2/(4.*MT**2))*cmath.log((1 + 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2))/(1 - 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2)))**2 + (-1 + MSd**2/(4.*MT**2))**2*cmath.log((1 + 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2))/(1 - 2*cmath.sqrt(((-1 + MSd**2/(4.*MT**2))*MT**2)/MSd**2)))**4))/MSd**8 )))/(cmath.pi**2*vev)',
                texname = 'G_S')

mu2h = Parameter(name = 'mu2h',
                 nature = 'internal',
                 type = 'real',
                 value = '-(lam1*vev**2) - (lam3*vevD**2)/2.',
                 texname = '\\mu _H')

mu2sd = Parameter(name = 'mu2sd',
                  nature = 'internal',
                  type = 'real',
                  value = '-0.5*(lam3*vev**2) - lam2*vevD**2',
                  texname = '\\mu _2')

I1e33 = Parameter(name = 'I1e33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I1e33}')

I2e33 = Parameter(name = 'I2e33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I2e33}')

I3e33 = Parameter(name = 'I3e33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I3e33}')

I4e33 = Parameter(name = 'I4e33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I4e33}')

