# Guide
# Pkg.dev("SixDOF")
# Pkg.add("SixDOF")
# Pkg.test("SixDOF")
# include("SixDOF.jl")

using Plots
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
m = 13.5
Ixx = .8244
Iyy = 1.135
Izz = 1.759
Ixz = 0.1204
mp = MassProp(m, Ixx, Iyy, Izz, Ixz)
nothing # hide
```

Next, we specify reference areas and lengths:


```
Sref = 0.55
bref = 2.8956
cref = 0.18994
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



CL0 = 0.28 # Zero-alpha lift
CLalpha = 3.45  # lift curve slope
CLq = 0.0 # Pitch rate derivative
CLM = 0.0 # Mach derivative
CLdf = 0.0  # flaps derivative
CLde = -0.36  # elevator derivative
CLmax = 1.4  # max CL (stall)
CLmin = -0.9  # min CL (neg stall)
alphas = 20*pi/180

CD0 = 0.03  # zero-lift drag coerff
U0 = 30.0  # velocity corresponding to Reynolds number of CD0
exp_Re = -0.2  # exponent in Reynolds number calc
e = 0.9  # Oswald efficiency
Mcc = 0.7  # crest critcal Mach number
CDdf = 0.0  # flaps
CDde = 0.0  # elevator
CDda = 0.0  # aileron
CDdr = 0.0  # rudder

CYbeta = -0.98 # Sideslip derivative
CYp = 0.0  # Roll rate derivative
CYr = 0.0 # Yaw rate derivative
CYda = 0.0 # Roll control (aileron) derivative
CYdr = -0.17 # Yaw control (rudder) derivative

Clbeta = -0.12  # Sideslip derivative
Clp = -0.26  # Roll rate derivative
Clr = 0.14  # Yaw rate derivative
Clda = 0.08  # Roll (aileron) control derivative
Cldr = 0.105  #Yaw (rudder) control derivative

Cm0 = -0.02338 # Zero-alpha pitch
Cmalpha = -0.38 # Alpha derivative
Cmq = -3.6 # Pitch rate derivative
CmM = 0.0
Cmdf = 0.0
Cmde = -0.5 # Pitch control derivative

Cnbeta = 0.25  # Slideslip derivative
Cnp = 0.022  # Roll rate derivative
Cnr = -0.35  # Yaw rate derivative
Cnda = 0.06  # Roll (aileron) control derivative
Cndr = 0.0  # Yaw (rudder) control derivative

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
alpha = 3.0*pi/180
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

leftWing = Wing([0.0 0.0 0.0],[-.1 -1 -.05], .3, .1, 0.0, 3.0)
rightWing = Wing([0.0 0.0 0.0],[-.1  1 -.05], .3, .1, 0.0, 3.0)

planeModel = (leftWing,rightWing)

animvtk(JACOBIAN, planeModel,"Aerosonde")
animvtk(NEWJACOBIAN, planeModel,"Aerosonde2")
