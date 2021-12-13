using VortexLattice

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [-1.0*pi/180, -1.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false

# horizontal stabilizer
xle_h = [0.0, 0.14]
yle_h = [0.0, 1.25]
zle_h = [0.0, 0.0]
chord_h = [0.7, 0.42]
theta_h = [0.0, 0.0]
phi_h = [0.0, 0.0]
fc_h = fill((xc) -> 0, 2) #camberline function for each section
ns_h = 6
nc_h = 3
spacing_s_h = Uniform()
spacing_c_h = Uniform()
mirror_h = false

# vertical stabilizer
xle_v = [0.0, 0.14]
yle_v = [0.0, 0.0]
zle_v = [0.0, 1.0]
chord_v = [0.7, 0.42]
theta_v = [0.0, 0.0]
phi_v = [0.0, 0.0]
fc_v = fill((xc) -> 0, 2) #camberline function for each section
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false

Sref = 8.0
cref = 0.9
bref = 10.0
rref = [1.9437, 0.0, 0.3708]#TODO Always make this the correct CoM
Vinf =  21.0967#Was originally 1
ref = Reference(Sref, cref, bref, rref, Vinf)

alpha = 1.2129*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = [false, false, false,false,false]

# generate surface panels for wing
wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)
wgrid2, wing2 = wing_to_surface_panels(xle, -yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

# generate surface panels for horizontal tail
hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(hgrid, [4.0, 0.0, 0.0])
translate!(htail, [4.0, 0.0, 0.0])
hgrid2, htail2 = wing_to_surface_panels(xle_h, -yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(hgrid2, [4.0, 0.0, 0.0])
translate!(htail2, [4.0, 0.0, 0.0])

# generate surface panels for vertical tail
vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vgrid, [4.0, 0.0, 0.0])
translate!(vtail, [4.0, 0.0, 0.0])

grids = [wgrid, hgrid,wgrid2, hgrid2, vgrid] #TODO
surfaces = [wing, htail,wing2, htail2, vtail] #TODO
surface_id = [1, 2, 3,4,5]

system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

CF, CM = body_forces(system; frame=Wind())

CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

dCF, dCM = stability_derivatives(system)

CDa, CYa, CLa = dCF.alpha
Cla, Cma, Cna = dCM.alpha
CDb, CYb, CLb = dCF.beta
Clb, Cmb, Cnb = dCM.beta
CDp, CYp, CLp = dCF.p
Clp, Cmp, Cnp = dCM.p
CDq, CYq, CLq = dCF.q
Clq, Cmq, Cnq = dCM.q
CDr, CYr, CLr = dCF.r
Clr, Cmr, Cnr = dCM.r


#Visualize the model
#properties= get_surface_properties(system)
#write_vtk("Anim/model",surfaces,properties; symmetric)

using SixDOF

CLM = 0.0 # Mach derivative
CLdf = 0.0  # flaps derivative
CLde = 0.0  # elevator derivative
alphas = 15*pi/180
U0=21.0967
exp_Re = -0.2  # exponent in Reynolds number calc
Mcc = 0.7  # crest critcal Mach number
e = 0.80  # Oswald efficiency#TODO Adjust
CDdf = 0.0  # flaps
CDde = 0.0  # elevator
CDda = 0.0  # aileron
CDdr = 0.0  # rudder
CYda = 0.0 # Roll control (aileron) derivative
CYdr = 0.0 # Yaw control (rudder) derivative
Clda = 0.0  # Roll (aileron) control derivative
Cldr = 0.0  #Yaw (rudder) control derivative
Cmdf = 0.0
Cmde = 0.0 # Pitch control derivative
CmM = 0.0
Cnda = 0.0  # Roll (aileron) control derivative
Cndr = 0.0  # Yaw (rudder) control derivative

m = 7.25755
Ixx = 38.05
Iyy = 26.46
Izz = 62.97
Ixz = -2.133
mp = MassProp(m, Ixx, Iyy, Izz, Ixz)

Wi = [0.0, 0.0, 0.0]
Wb = [0.0, 0.0, 0.0]
rho = 1.225
asound = 300.0
g = 9.81
atm = ConstantAtmosphere(Wi, Wb, rho, asound, g)

controller = ConstantController(0.0, 0.0, 0.0, 0.0, 0.0)

sd = StabilityDeriv(CL, CLa, CLq, CLM, CLdf, CLde, alphas,
    CD, U0, exp_Re, e, Mcc, CDdf, CDde, CDda, CDdr,
    CYb, -CYp, -CYr, CYda, CYdr,
    -Clb, Clp, Clr, Clda, Cldr,
    Cm, Cma, Cmq, CmM, Cmdf, Cmde,
    -Cnb, Cnp, Cnr, Cnda, Cndr) #some are negative to account for changes in refrence frame.
nothing # hide

CT0 = 0.11221
CT1 = -0.13803
CT2 = -0.047394
CQ0 = 0.0062
CQ1 = 0.00314
CQ2 = -0.015729
D = 10*0.0254
num = 0
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
alpha = 1.21290*pi/180
s0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, U0, 0.0, 0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 4.0)
p = mp, ref, sd, propulsion, inertial2, atm, controller


#Tyson testing his stuff
JACOBIAN=getjacobian(s0,p)
s0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, U0, 0.0, 0.0, 0.0, 0.0, 0.0]
p = mp, ref, sd, propulsion, inertial2, atm, controller

leftWing = Wing([0.0 0.0 0.0],[-.1 -1 -.2], .3, .1, 0.0, 3.0)
rightWing = Wing([0.0 0.0 0.0],[-.1  1 -.2], .3, .1, 0.0, 3.0)
leftTail = Wing([-0.8 0.0 0.0],[-.05 -.33 -0.0], .3, .05, 0.0, 0.0)
rightTail = Wing([-0.8 0.0 0.0],[-.05  .33 -0.0], .3, .05, 0.0, 0.0)
rudder = Wing([-0.8 0.0 0.0],[-.05  0.0 -0.33], .3, .05, 0.0, 0.0)
planeModel = (leftWing,rightWing,leftTail,rightTail,rudder)
animvtk(JACOBIAN, planeModel,"Test")


using LinearAlgebra
Eigenvals, Eigenvecs =eigen(JACOBIAN)

longitudinal = [JACOBIAN[7,7] JACOBIAN[7,9] JACOBIAN[7,11] JACOBIAN[7,5]
                JACOBIAN[9,7] JACOBIAN[9,9] JACOBIAN[9,11] JACOBIAN[9,5]
                JACOBIAN[11,7] JACOBIAN[11,9] JACOBIAN[11,11] JACOBIAN[11,5]
                JACOBIAN[5,7] JACOBIAN[5,9] JACOBIAN[5,11] JACOBIAN[5,5]]
EigenvalsLo, EigenvecsLo =eigen(longitudinal)
lateral = [JACOBIAN[8,8] JACOBIAN[8,10] JACOBIAN[8,12] JACOBIAN[8,4]
                JACOBIAN[10,8] JACOBIAN[10,10] JACOBIAN[10,12] JACOBIAN[10,4]
                JACOBIAN[12,8] JACOBIAN[12,10] JACOBIAN[12,12] JACOBIAN[12,4]
                JACOBIAN[4,8] JACOBIAN[4,10] JACOBIAN[4,12] JACOBIAN[4,4]]
EigenvalsLa, EigenvecsLa =eigen(lateral)


#Hunting down errors--------------------------------
s1=s0
state = State(s0...)
Fa,Ma = aeroforces(sd, atm, state, controller, ref, mp)
Fp,Mp = propulsionforces(propulsion, atm, state, controller, ref, mp)
Fg,Mg = gravityforces(inertial2, atm, state, controller, ref, mp)
ds = zeros(length(s0))
ds[1:3]=s0[7:9]
#required for use of sixdof!() Not implemented.
time = 0.0
#Also must be reset each time sixdof!() is used.
state0=s0
params0=p
sixdof!(ds, state0, params0, time)#sixdof!() modifies ds so that it is now the defailt derivatives
#Result: values in ds give the derivatives about the center point.
dsc=ds

s1[10]=s1[10]+.001
state1 = State(s1...)
Fap,Map = aeroforces(sd, atm, state1, controller, ref, mp)
Fpp,Mpp = propulsionforces(propulsion, atm, state1, controller, ref, mp)
Fgp,Mgp = gravityforces(inertial2, atm, state1, controller, ref, mp)
#F = Fa + Fp + Fg
time = 0.0
state0=s0
params0=p
ds = zeros(length(s0))
ds[1:3]=s0[7:9]
#Perturb s at correct index
state0[10]=state0[10]+.001
#Run sixdof!() for perturbed value.
sixdof!(ds, state0, params0, time)
dsp=ds#perturbed ds

col11=(dsp-dsc)/.001
