# This file was automatically created by FeynRules 2.3.49
# Mathematica version: 12.1.0 for Linux x86 (64-bit) (March 18, 2020)
# Date: Wed 23 Oct 2024 21:15:50



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
sina = Parameter(name = 'sina',
                 nature = 'external',
                 type = 'real',
                 value = 0.2,
                 texname = '\\text{sina}',
                 lhablock = 'NPINPUTS',
                 lhacode = [ 1 ])

MQloop = Parameter(name = 'MQloop',
                   nature = 'external',
                   type = 'complex',
                   value = MT,
                   texname = '\\text{MQloop}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 2 ])

AAC = Parameter(name = 'AAC',
                nature = 'external',
                type = 'real',
                value = 1.e-6,
                texname = 'G_{\\text{aachi}}',
                lhablock = 'NPINPUTS',
                lhacode = [ 3 ])

ychi1 = Parameter(name = 'ychi1',
                  nature = 'external',
                  type = 'real',
                  value = 1.,
                  texname = 'y_{\\text{chi1}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 4 ])

ychi0 = Parameter(name = 'ychi0',
                  nature = 'external',
                  type = 'real',
                  value = 0.,
                  texname = 'y_{\\text{chi0}}',
                  lhablock = 'NPINPUTS',
                  lhacode = [ 5 ])

ychi10 = Parameter(name = 'ychi10',
                   nature = 'external',
                   type = 'real',
                   value = 0.,
                   texname = 'y_{\\text{chi10}}',
                   lhablock = 'NPINPUTS',
                   lhacode = [ 6 ])

ystt = Parameter(name = 'ystt',
                 nature = 'external',
                 type = 'real',
                 value = (0.28284271247461906*MT)/vev,
                 texname = 'y_{\\text{Sdtt}}',
                 lhablock = 'NPINPUTS',
                 lhacode = [ 7 ])

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

M1 = Parameter(name = 'M1',
               nature = 'external',
               type = 'real',
               value = 500.,
               texname = '\\text{M1}',
               lhablock = 'MASS',
               lhacode = [ 5000023 ])

M0 = Parameter(name = 'M0',
               nature = 'external',
               type = 'real',
               value = 425.,
               texname = '\\text{M0}',
               lhablock = 'MASS',
               lhacode = [ 5000022 ])

MSd = Parameter(name = 'MSd',
                nature = 'external',
                type = 'real',
                value = 1000,
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

w1 = Parameter(name = 'w1',
               nature = 'external',
               type = 'real',
               value = 0.,
               texname = '\\text{w1}',
               lhablock = 'DECAY',
               lhacode = [ 5000023 ])

w0 = Parameter(name = 'w0',
               nature = 'external',
               type = 'real',
               value = 0.,
               texname = '\\text{w0}',
               lhablock = 'DECAY',
               lhacode = [ 5000022 ])

WS = Parameter(name = 'WS',
               nature = 'external',
               type = 'real',
               value = 1.,
               texname = '\\text{WS}',
               lhablock = 'DECAY',
               lhacode = [ 55 ])

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

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

GGS = Parameter(name = 'GGS',
                nature = 'internal',
                type = 'complex',
                value = '-(( 1 + (7*MSd**2)/(120.*MQloop**2) + MSd**4/(168.*MQloop**4) + (13*MSd**6)/(16800.*MQloop**6) + (2*MSd**8)/(17325.*MQloop**8) + (19*MSd**10)/(1.009008e6*MQloop**10) + MSd**12/(305760.*MQloop**12) + (5*MSd**14)/(8.401536e6*MQloop**14) + (7*MSd**16)/(6.235515e7*MQloop**16) if MSd**2/(4.*MQloop**2)<1 else (1136*MQloop**16)/(3.*MSd**16) + (4003*MQloop**14)/(30.*MSd**14) + (49*MQloop**12)/MSd**12 + (37*MQloop**10)/(2.*MSd**10) + (6*MQloop**8)/MSd**8 - (6*MQloop**6)/MSd**6 + (6*MQloop**2)/MSd**2 - (2640*MQloop**16*cmath.log(-(MSd**2/MQloop**2)))/(7.*MSd**16) - (714*MQloop**14*cmath.log(-(MSd**2/MQloop**2)))/(5.*MSd**14) - (294*MQloop**12*cmath.log(-(MSd**2/MQloop**2)))/(5.*MSd**12) - (55*MQloop**10*cmath.log(-(MSd**2/MQloop**2)))/(2.*MSd**10) - (16*MQloop**8*cmath.log(-(MSd**2/MQloop**2)))/MSd**8 - (15*MQloop**6*cmath.log(-(MSd**2/MQloop**2)))/MSd**6 + (6*MQloop**4*cmath.log(-(MSd**2/MQloop**2)))/MSd**4 + (6*MQloop**4*cmath.log(-(MSd**2/MQloop**2))**2)/MSd**4 - (3*MQloop**2*cmath.log(-(MSd**2/MQloop**2))**2)/(2.*MSd**2) )*G**2*ystt)/(12.*cmath.pi**2*MQloop*cmath.sqrt(2))',
                texname = 'G_S')

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
                value = '-(G**2*(1 + (13*MH**6)/(16800.*MT**6) + MH**4/(168.*MT**4) + (7*MH**2)/(120.*MT**2)))/(12.*cmath.pi**2*vev)',
                texname = 'G_H')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

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

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

I1a33 = Parameter(name = 'I1a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I1a33}')

I2a33 = Parameter(name = 'I2a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I2a33}')

I3a33 = Parameter(name = 'I3a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I3a33}')

I4a33 = Parameter(name = 'I4a33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I4a33}')

