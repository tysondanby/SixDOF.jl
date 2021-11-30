
```

First, we import the module.

```
using SixDOF
```

There are several structs we need to define.  The inputs used in this example correspond to the Zagi flying wing in Appendix E of Small Unmanned Aircraft: Theory and Practice by Beard and McLain.  We specify the mass properties:

```
```

For this example:
```
m = 6.804
Ixx = 0.54711
Iyy = 2.5767834
Izz = 3.03369
Ixz = 0.240388
mp = MassProp(m, Ixx, Iyy, Izz, Ixz)
nothing # hide
```

Next, we specify reference areas and lengths:


```
Sref = 2.667
bref = 7.111
cref = 0.41666
ref = Reference(Sref, bref, cref)
nothing # hide

controller = ConstantController(0.0, 0.0, 0.0, 0.0, 0.8)
nothing # hide

Wi = [0.0, 0.0, 0.0]
Wb = [0.0, 0.0, 0.0]
rho = 1.2682
asound = 300.0
g = 9.81
atm = ConstantAtmosphere(Wi, Wb, rho, asound, g)
nothing # hide



CL0 = 0.10889 # Zero-alpha lift
CLalpha = 5.37212  # lift curve slope
CLq = 9.60847 # Pitch rate derivative INV
CLM = 0.0 # Mach derivative
CLdf = 0.0  # flaps derivative
CLde = 0.0  # elevator derivative
CLmax = 1.4  # max CL (stall)TODO
CLmin = -0.9  # min CL (neg stall)TODO
alphas = 20*pi/180 #TODO

CD0 = 0.00966  # zero-lift drag coerff
U0 = 34.68  # velocity corresponding to Reynolds number of CD0 TODO
exp_Re = -0.2  # exponent in Reynolds number calc TODO
e = 0.8  # Oswald efficiency TODO
Mcc = 0.7  # crest critcal Mach number TODO
CDdf = 0.0  # flaps TODO
CDde = 0.0  # elevator TODO
CDda = 0.0  # aileron TODO
CDdr = 0.0  # rudder TODO

CYbeta = -0.39446 # Sideslip derivative
CYp = .06773  # Roll rate derivative inv
CYr = -0.4642 # Yaw rate derivative inv
CYda = 0.0 # Roll control (aileron) derivative
CYdr = 0.0 # Yaw control (rudder) derivative

Clbeta = 0.03442  # Sideslip derivative  inv
Clp = -0.42453  # Roll rate derivative
Clr = 0.06703  # Yaw rate derivative
Clda = 0.0  # Roll (aileron) control derivative TODO
Cldr = 0.0  #Yaw (rudder) control derivative

Cm0 = -0.00212 # Zero-alpha pitch
Cmalpha = -1.99091 # Alpha derivative
Cmq = -24.97412 # Pitch rate derivative
CmM = 0.0 #TODO
Cmdf = 0.0 #TODO
Cmde = -0.3254 # Pitch control derivative TODO

Cnbeta = 0.21874  # Slideslip derivative Perhaps is inv?
Cnp = 0.02522  # Roll rate derivative
Cnr = -0.25595  # Yaw rate derivative
Cnda = -0.00328  # Roll (aileron) control derivative TODO
Cndr = 0.0  # Yaw (rudder) control derivative TODO

sd = StabilityDeriv(CL0, CLalpha, CLq, CLM, CLdf, CLde, alphas,
    CD0, U0, exp_Re, e, Mcc, CDdf, CDde, CDda, CDdr,
    CYbeta, CYp, CYr, CYda, CYdr,
    Clbeta, Clp, Clr, Clda, Cldr,
    Cm0, Cmalpha, Cmq, CmM, Cmdf, Cmde,
    Cnbeta, Cnp, Cnr, Cnda, Cndr)
nothing # hide

CT0 = 0.11221
CT1 = -0.13803
CT2 = -0.047394
CQ0 = 0.0062
CQ1 = 0.00314
CQ2 = -0.015729
D = 10*0.0254
num = 2
type = COCOUNTER
R = 0.5
Kv = 2500.0 * pi/30
i0 = 0.3
voltage = 8.0
propulsion = MotorPropBatteryDataFit(CT2, CT1, CT0, CQ2, CQ1, CQ0, D, num, type, R, Kv, i0, voltage)
nothing # hide

inertial2 = UniformGravitationalField()
nothing # hide


Vinf = U0
alpha = -0.52*pi/180
s0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Vinf*cos(alpha), 0.0, Vinf*sin(alpha), 0.0, 0.0, 0.0]
tspan = (0.0, 4.0)
p = mp, ref, sd, propulsion, inertial2, atm, controller

# using DifferentialEquations
# prob = DifferentialEquations.ODEProblem(sixdof!, s0, tspan, p)
# sol = DifferentialEquations.solve(prob)
# nothing # hide

#Tyson testing his stuff
JACOBIAN=finite_jacobian(s0,p)
s0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Vinf*cos(alpha), 0.0, Vinf*sin(alpha), 0.0, 0.0, 0.0]
p = mp, ref, sd, propulsion, inertial2, atm, controller
NEWJACOBIAN=getjacobian(s0,p)

#=
struct wing{TV, TF}
    pos::TV
    dir::TV
    cr::TF
    ct::TF
    tilt::TF
    twist::TF
end
=#

leftWing = Wing([0.0 0.0 0.0],[-.1 -1.5 0.0], .4, .1, 0.0, 0.0)
rightWing = Wing([0.0 0.0 0.0],[-.1  1.5 0.0], .4, .1, 0.0, 0.0)
leftTail = Wing([-2.5 0.0 0.0],[-.03 -.75 0.0], .3, .1, 0.0, 0.0)
rightTail = Wing([-2.5 0.0 0.0],[-.03  .75 0.0], .3, .1, 0.0, 0.0)
rudder = Wing([-2.5 0.0 0.0],[-0.1  0.0 -0.75], .3, .1, 0.0, 0.0)

planeModel = (leftWing,rightWing,leftTail,rightTail,rudder)

animvtk(JACOBIAN, planeModel,"DBF")
animvtk(NEWJACOBIAN, planeModel,"DBF2")
