
module SixDOF


using ForwardDiff
using LinearAlgebra: norm, cross, eigen
using WriteVTK


export Control, MassProp, Reference
export AbstractAeroModel, AbstractPropulsionModel, AbstractInertialModel,
    AbstractAtmosphereModel, AbstractController
export StabilityDeriv, MotorPropBatteryDataFit, UniformGravitationalField,
    ConstantAtmosphere, ConstantController
export CO, COUNTER, COCOUNTER
export sixdof!
export finite_jacobian, geteig, getjacobian, animvtk
export Wing
#DEBUG
export aeroforces,propulsionforces,gravityforces, State


# ------ General Structs -------

struct Wing{TV, TF}
    pos::TV
    dir::TV
    cr::TF
    ct::TF
    tilt::TF
    twist::TF
end

"""
    State(x, y, z, phi, theta, psi, u, v, w, p, q, r)

State of the aircraft: positions in inertial frame, euler angles,
velocities in body frame, angular velocities in body frame.
"""
struct State{TF}
    x::TF  # position (inertial frame)
    y::TF
    z::TF
    phi::TF  # orientation, euler angles
    theta::TF
    psi::TF
    u::TF  # velocity (body frame)
    v::TF
    w::TF
    p::TF  # angular velocity (body frame)
    q::TF
    r::TF
end

"""
    Control(de, dr, da, df, throttle)

Define the control settings: delta elevator, delta rudder, delta aileron,
delta flaps, and throttle.
"""
struct Control{TF}
    de::TF  # elevator
    dr::TF  # rudder
    da::TF  # aileron
    df::TF  # rudder
    throttle::TF
end

"""
    MassProp(m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz)

Mass and moments of inertia in the body frame.
Ixx = int(y^2 + z^2, dm)
Ixz = int(xz, dm)

Most aircraft are symmetric about y and so there is a convenience method
to specify only the nonzero components.
MassProp(m, Ixx, Iyy, Izz, Ixz)
"""
struct MassProp{TF}
    m::TF
    Ixx::TF
    Iyy::TF
    Izz::TF
    Ixz::TF
    Ixy::TF
    Iyz::TF
end

# most aircraft are symmetric in y
MassProp(m, Ixx, Iyy, Izz, Ixz) = MassProp(m, Ixx, Iyy, Izz, Ixz, zero(Ixx), zero(Ixx))


"""
    Reference(S, b, c)

The reference area, span, and chord used in the aerodynamic computations.
"""
struct Reference{TF}
    S::TF  # area
    b::TF  # span
    c::TF  # chord
end

# ----------------------------------------------


# --------------- Interfaces ---------------

abstract type AbstractAtmosphereModel end

"""
    wind(model::AbstractAtmosphereModel, state)

Compute wind velocities.

**Returns**
- Wi: wind velocities in inertial frame
- Wb: gust velocities in body frame (just a convenience to allow some velocities in body frame)
"""
function wind(model::AbstractAtmosphereModel, state)
    @warn "wind function not implemented for AbstractAtmosphereModel"
    Wi = [0.0, 0.0, 0.0]
    Wb = [0.0, 0.0, 0.0]
    return Wi, Wb
end

"""
    properties(model::AbstractAtmosphereModel, state)

Compute atmospheric density and the speed of sound.
"""
function properties(model::AbstractAtmosphereModel, state)
    @warn "properties function not implemented for AbstractAtmosphereModel"
    rho = 1.225  # sea-level properties
    asound = 340.3
    return rho, asound
end

"""
    gravity(model::AbstractAtmosphereModel, state)

Compute the local acceleration of gravity.
"""
function gravity(model::AbstractAtmosphereModel, state)
    @warn "gravity function not implemented for AbstractAtmosphereModel"
    g = 9.81
    return g
end


# ----

abstract type AbstractAeroModel end

"""
    aeroforces(model::AbstractAeroModel, atm::AbstractAtmosphereModel, state::State, control::Control, mp::MassProp, ref::Reference)

Compute the aerodynamic forces and moments in the body frame.
return F, M
"""
function aeroforces(model::AbstractAeroModel, atm, state, control, mp, ref)
    @warn "aeroforces function not implemented for AbstractAeroModel"
    # forces and moments in body frame
    F = [0.0, 0.0, 0.0]
    M = [0.0, 0.0, 0.0]
    return F, M
end

# ----

abstract type AbstractPropulsionModel end

"""
    propulsionforces(model::AbstractPropulsionModel, atm::AbstractAtmosphereModel, state::State, control::Control, mp::MassProp, ref::Reference)

Compute the propulsive forces and moments in the body frame.
return F, M
"""
function propulsionforces(model::AbstractPropulsionModel, atm, state, control, mp, ref)
    @warn "propulsionforces function not implemented for AbstractPropulsionModel"
    # forces and moments in body frame
    F = [0.0, 0.0, 0.0]
    M = [0.0, 0.0, 0.0]
    return F, M
end

# ----

abstract type AbstractInertialModel end

"""
    gravityforces(model::AbstractInertialModel, atm::AbstractAtmosphereModel, state::State, control::Control, mp::MassProp, ref::Reference)

Compute the gravitational forces and moments in the body frame.
return F, M
"""
function gravityforces(model::AbstractInertialModel, atm, state, control, mp, ref)
    @warn "gravityforces function not implemented for AbstractInertialModel"
    # forces and moments in body frame
    F = [0.0, 0.0, 0.0]
    M = [0.0, 0.0, 0.0]
    return F, M
end


# ----

abstract type AbstractController end

"""
    setcontrol(controller::AbstractController, time, atm::AbstractAtmosphereModel, state::State, lastcontrol::Control, mp::MassProp, ref::Reference)

Compute control state for next time step given current state
return control::Control
"""
function setcontrol(controller::AbstractController, time, atm, state, mp, ref)
    @warn "setcontrol function not implemented for AbstractController"
    control = Control(0.0, 0.0, 0.0, 0.0, 0.0)
    return control
end


# -----------------------------




# ------------- helper functions (private) --------------


"""
    inertialtobody(state)

Construct a rotation matrix from inertial frame to body frame

The assumed order of rotation is
1) psi radians about the z axis,
2) theta radians about the y axis,
3) phi radians about the x axis.

This is an orthogonal transformation so its inverse is its transpose.
"""
function inertialtobody(state)

    R = Array{eltype(state.phi)}(undef, 3, 3)

    cphi, ctht, cpsi = cos.([state.phi, state.theta, state.psi])
    sphi, stht, spsi = sin.([state.phi, state.theta, state.psi])

    R[1, 1] = ctht*cpsi
    R[1, 2] = ctht*spsi
    R[1, 3] = -stht

    R[2, 1] = sphi*stht*cpsi - cphi*spsi
    R[2, 2] = sphi*stht*spsi + cphi*cpsi
    R[2, 3] = sphi*ctht

    R[3, 1] = cphi*stht*cpsi + sphi*spsi
    R[3, 2] = cphi*stht*spsi - sphi*cpsi
    R[3, 3] = cphi*ctht

    return R

end


"""
    windtobody(alpha, beta)

Rotation matrix from wind frame to body frame.
- alpha: angle of attack
- beta: sideslip angle
"""
function windtobody(alpha, beta)

    ca, cb = cos.([alpha, beta])
    sa, sb = sin.([alpha, beta])

    Rwb = [ca*cb  -ca*sb  -sa;
           sb      cb     0.0;
           sa*cb  -sa*sb  ca]

    return Rwb
end

"""
    stabilitytobody(alpha, beta)

Rotation matrix from stability frame to body frame.
- alpha: angle of attack
"""
function stabilitytobody(alpha)

    ca = cos(alpha)
    sa = sin(alpha)

    Rsb = [ca  0.0  -sa;
           0.0 1.0  0.0;
           sa  0.0  ca]

    return Rsb
end


"""
    windaxes(atm::AbstractAtmosphereModel, state)

Compute relative velocity in wind axes (airspeed, aoa, sideslip)
"""
function windaxes(atm, state)

    # velocity vectors
    Vb = [state.u, state.v, state.w]
    Wi, Wb = wind(atm, state)

    Rib = inertialtobody(state)

    # relative wind
    Vrel = Vb - (Rib*Wi + Wb)

    # airspeed
    Va = norm(Vrel)

    # angle of attack
    alpha = atan(Vrel[3], Vrel[1])

    # sideslip
    beta = asin(Vrel[2] / Va)

    return Va, alpha, beta

end


# ----------------------------------------------------


# ------- Some Default Interface Implementations -----

"""
    StabilityDeriv(CL0, CLalpha, CLq, CLM, CLdf, CLde, alphas,
        CD0, U0, exp_Re, e, Mcc, CDdf, CDde, CDda, CDdr,
        CYbeta, CYp, CYr, CYda, CYdr, Clbeta,
        Clp, Clr, Clda, Cldr,
        Cm0, Cmalpha, Cmq, CmM, Cmdf, Cmde,
        Cnbeta, Cnp, Cnr, Cnda, Cndr)

Stability derivatives of the aircraft.  Most are self explanatory if you are
familiar with stability derivatives (e.g., CLalpha is dCL/dalpha or the lift curve slope).
Some less familiar ones include
- M: Mach number
- alphas: the angle of attack for stall
- U0: the speed for the reference Reynolds number CD0 was computed at
- exp_Re: the coefficient in the denominator of the skin friction coefficient (0.5 laminar, 0.2 turbulent)
- e: Oswald efficiency factor
- Mcc: crest critical Mach number (when compressibility drag rise starts)

"""
struct StabilityDeriv{TF} <: AbstractAeroModel
    CL0::TF
    CLalpha::TF
    CLq::TF
    CLM::TF
    CLdf::TF
    CLde::TF
    alphas::TF  # TODO: should probably do in terms of CLmax

    CD0::TF
    U0::TF  # velocity corresponding to Reynolds number of CD0  (TODO: rethink this)
    exp_Re::TF  # exponent for Reynolds number scaling. typical values: exp_Re = 0.5 laminar, 0.2 turbulent
    e::TF  # Oswald efficiency factor
    Mcc::TF  # crest-critical Mach number when compressibility drag rise starts (quartic)
    CDdf::TF
    CDde::TF
    CDda::TF
    CDdr::TF

    CYbeta::TF
    CYp::TF
    CYr::TF
    CYda::TF
    CYdr::TF

    Clbeta::TF
    Clp::TF
    Clr::TF
    Clda::TF
    Cldr::TF

    Cm0::TF
    Cmalpha::TF
    Cmq::TF
    CmM::TF
    Cmdf::TF
    Cmde::TF

    Cnbeta::TF
    Cnp::TF
    Cnr::TF
    Cnda::TF
    Cndr::TF
end


"""
A simple (mostly) linear aerodynamics model
"""
function aeroforces(sd::StabilityDeriv, atm, state, control, ref, mp)

    # airspeed, angle of attack, sideslip
    Va, alpha, beta = windaxes(atm, state)

    # Mach number and dynamic pressure
    rho, asound = properties(atm, state)
    Mach = Va / asound
    qdyn = 0.5 * rho * Va^2

    # rename for convenience
    p = state.p
    q = state.q
    r = state.r
    de = control.de
    df = control.df
    dr = control.dr
    da = control.da


    # lift
    CL = sd.CL0 + sd.CLalpha*alpha + sd.CLq*q *ref.c/(2*Va) + sd.CLM*Mach
        + sd.CLdf*df + sd.CLde*de

    em = exp(-50*(alpha - sd.alphas))
    ep = exp(50*(alpha + sd.alphas))
    sigma = (1 + em + ep)/((1 + em)*(1 + ep))
    CL = (1- sigma)*CL + sigma * 2 * sign(alpha)*sin(alpha)^2*cos(alpha)

    # drag
    CDp = sd.CD0*(Va/sd.U0)^sd.exp_Re
    CDi = CL^2/(pi*(ref.b^2/ref.S)*sd.e)
    CDc = (Mach < sd.Mcc) ? 0.0 : 20*(Mach - sd.Mcc)^4

    CD = CDp + CDi + CDc + abs(sd.CDdf*df) + abs(sd.CDde*de) + abs(sd.CDda*da) + abs(sd.CDdr*dr)

    # side force
    CY = sd.CYbeta*beta + (sd.CYp*p + sd.CYr*r)*ref.b/(2*Va) + sd.CYda*da + sd.CYdr*dr

    # rolling moment
    Cl = sd.Clbeta*beta + (sd.Clp*p + sd.Clr*r)*ref.b/(2*Va) + sd.Clda*da + sd.Cldr*dr

    # pitching moment
    Cm = sd.Cm0 + sd.Cmalpha*alpha + sd.Cmq*q * ref.c/(2*Va) + sd.CmM*Mach + sd.Cmdf*df + sd.Cmde*de

    # yawing moment
    Cn = sd.Cnbeta*beta + (sd.Cnp*p + sd.Cnr*r)*ref.b/(2*Va) + sd.Cnda*da + sd.Cndr*dr

    # transfer forces from wind to body axes
    Rwb = windtobody(alpha, beta)

    F = Rwb*[-CD, CY, -CL] * qdyn * ref.S

    M = Rwb*[Cl*ref.b, Cm*ref.c, Cn*ref.b] * qdyn * ref.S

    return F, M
end

@enum PropType CO=1 COUNTER=-1 COCOUNTER=0

"""
    MotorPropBatteryDataFit(CT2, CT1, CT0, CQ2, CQ1, CQ0, D, num, type,
        R, Kv, i0, voltage)

**Inputs**
- CT2, CT1, CT0: quadratic fit to propeller thrust coefficient of form: CT = CT2*J2 + CT1*J + CT0
- CQ2, CQ1, CQ0: quadratic fit to propeller torque coefficient of form: CQ = CQ2*J2 + CQ1*J + CQ0
- D: propeller diameter
- num: number of propellers
- type: CO (torques add), COUNTER (torques add but with minus sign), COCOUNTER (no torque, they cancel out)
- R: motor resistance
- Kv: motor Kv
- i0: motor no-load current
- voltage: battery voltage
"""
struct MotorPropBatteryDataFit{TF, TI, PropType} <: AbstractPropulsionModel
    # CT = CT2*J2 + CT1*J + CT0
    # CQ = CQ2*J2 + CQ1*J + CQ0
    CT2::TF  # prop data fit
    CT1::TF
    CT0::TF
    CQ2::TF
    CQ1::TF
    CQ0::TF
    D::TF  # prop diameter
    num::TI
    type::PropType

    R::TF  # motor resistance
    Kv::TF  # motor Kv
    i0::TF  # motor no-load current

    voltage::TF  # battery voltage
end

function propulsionforces(prop::MotorPropBatteryDataFit, atm, state, control, ref, mp)

    # airspeed, angle of attack, sideslip
    Va, _, _ = windaxes(atm, state)

    # density
    rho, _ = properties(atm, state)

    D = prop.D

    # determine torque for motor/prop match (quadratic equation)
    a = rho*D^5/(2*pi)^2 * prop.CQ0
    b = rho*D^4/(2*pi)*Va * prop.CQ1 + 1.0/(prop.R*prop.Kv)
    c = rho*D^3*Va^2 * prop.CQ2 - control.throttle*prop.voltage/(prop.R*prop.Kv) + prop.i0/prop.Kv
    Omega = (-b + sqrt(b^2 - 4*a*c))/(2*a)

    # advance ratio
    n = Omega/(2*pi)
    J = Va/(n*D)

    # thrust and torque
    CT = prop.CT0 + prop.CT1*J + prop.CT2*J^2
    CQ = prop.CQ0 + prop.CQ1*J + prop.CQ2*J^2

    T = prop.num * CT * rho * n^2 * D^4
    Q = prop.num * CQ * rho * n^2 * D^5 * Int(prop.type)

    return [T, 0, 0], [Q, 0, 0]
end

"""
    UniformGravitationalField()

Assumes center of mass and center of gravity are coincident.
"""
struct UniformGravitationalField <: AbstractInertialModel end

function gravityforces(model::UniformGravitationalField, atm, state, control, ref, mp)

    W = mp.m * gravity(atm, state)
    ct, cp = cos.([state.theta, state.phi])
    st, sp = sin.([state.theta, state.phi])

    Fg = W*[-st, ct*sp, ct*cp]
    Mg = [zero(W), zero(W), zero(W)]  # no gravitational moment

    return Fg, Mg
end


"""
    ConstantAtmosphere(Wi, Wb, rho, asound, g)

Constant atmospheric properties.
"""
struct ConstantAtmosphere{TF, TV<:AbstractVector{TF}} <: AbstractAtmosphereModel
    Wi::TV
    Wb::TV
    rho::TF
    asound::TF
    g::TF
end


function wind(atm::ConstantAtmosphere, state)
    return atm.Wi, atm.Wb
end

function properties(atm::ConstantAtmosphere, state)
    return atm.rho, atm.asound
end

function gravity(atm::ConstantAtmosphere, state)
    return atm.g
end


"""
    ConstantController(de, dr, da, df, throttle)

Just a dummy controller that outputs constant control outputs the whole time.
"""
struct ConstantController{TF} <: AbstractController
    de::TF
    dr::TF
    da::TF
    df::TF
    throttle::TF
end

function setcontrol(controller::ConstantController, time, atm, state, mp, ref)
    return Control(controller.de, controller.dr, controller.da, controller.df, controller.throttle)
end
# --------------------------------------------------------


# ------------- main functions (public) --------------

"""
    sixdof!(ds, s, params, time)

dynamic and kinematic ODEs.  Follows format used in DifferentialEquations package.
- s = x, y, z, phi, theta, psi, u, v, w, p, q, r (same order as State)
- params = control, massproperties, reference, aeromodel, propmodel, inertialmodel, atmmodel
"""
function sixdof!(ds, s, params, time)

    x, y, z, phi, theta, psi, u, v, w, p, q, r = s
    mp, ref, aeromodel, propmodel, inertialmodel, atmmodel, controller = params

    # ---- controller -------
    state = State(s...)
    control = setcontrol(controller, time, atmmodel, state, mp, ref)
    # -----------------------

    # --------- forces and moments ---------
    # aerodynamics
    Fa, Ma = aeroforces(aeromodel, atmmodel, state, control, ref, mp)

    # propulsion
    Fp, Mp = propulsionforces(propmodel, atmmodel, state, control, ref, mp)

    # weight
    Fg, Mg = gravityforces(inertialmodel, atmmodel, state, control, ref, mp)

    # total forces and moments
    F = Fa + Fp + Fg
    M = Ma + Mp + Mg

    # --------------------------------------


    # ----- derivative of state --------
    Vb = [u, v, w]
    omegab = [p, q, r]

    # linear kinematics
    Rib = inertialtobody(state)
    rdot = Rib' * Vb

    # angular kinematics
    phidot = p + (q*sin(phi) + r*cos(phi))*tan(theta)
    thetadot = q*cos(phi) - r*sin(phi)
    psidot = (q*sin(phi) + r*cos(phi))/cos(theta)

    # linear dynamics
    vdot = F/mp.m - cross(omegab, Vb)

    # angular dyna1mics
    I = [mp.Ixx -mp.Ixy -mp.Ixz;
         -mp.Iyz mp.Iyy -mp.Iyz;
         -mp.Ixz -mp.Iyz mp.Izz]
    omegadot = I \ (M - cross(omegab, I*omegab))

    # -------------------------

    # TODO: if we need more efficiency we can avoid allocating and then assigning.
    ds[1:3] = rdot
    #ds[4:6] = omegab #Use for quaternions??
    ds[4] = phidot
    ds[5] = thetadot
    ds[6] = psidot
    ds[7:9] = vdot
    ds[10:12] = omegadot
end


function multQuat(a,b)#a and b are quaternions (four element vectors) being multiplied such that
    #result = a*b
    w = a[1]*b[1] - b[2]*a[2] - b[3]*a[3] - b[4]*a[4]
    x = b[1]*a[2] + b[2]*a[1] - b[3]*a[4] + b[4]*a[3]
    y = b[1]*a[3] + b[2]*a[4] + b[3]*a[1] - b[4]*a[2]
    z = b[1]*a[4] - b[2]*a[3] + b[3]*a[2] + b[4]*a[1]
    out = [w,x,y,z]
    return out
end

#Rotates a quaternion (quat) that represents a rigid body.
#Rotation is based on the angular velocity (omeg) and the time step (dt).
function rotQuat(p,omeg,dt) #(4-tuple,3-tuple,double)
    #result = a rotated quaternion representing the rigid body one dt later.

    #for use later
    oriMag= norm(p)

    #find axis of rotation and angle to rotate through.
    angSpd= norm(omeg)#angular speed about dir

    if angSpd == 0.0 #AVOID NaNs
        dir= omeg
    else
        dir= omeg ./ angSpd #Axis of rotation
    end

    th= angSpd*dt#theta (angular displacement) about the dir axis

    #create quaternions for rotation
    th = th/2.0 #quaternions double the angle input to them
    q= [  cos(th), sin(th)*dir[1], sin(th)*dir[2], sin(th)*dir[3]  ]#the rotation quaternion
    qinv= [  cos(-th), sin(-th)*dir[1], sin(-th)*dir[2], sin(-th)*dir[3]  ]#other rotation quaternion

    #make output
    out= multQuat(multQuat(q,p),qinv)
    finalMag= norm(out) #Make sure magnitude is 1

    if finalMag != 0
        out = out.*(oriMag/finalMag)#DONT NORMALIZE THE QUATERNION!! Make the output have the same magnitude as before.
    end

    return out
end

function eulerToQuat(phi,th,psi)#takes euler angles. Gives axis and angle.
    #Define quaternions for the whole rotation.
    phi = phi/2.0
    th = th/2.0
    psi = psi/2.0

    q1= [cos(psi),0,0,sin(psi)]
    yb= rotQuat([0.0,0.0,1.0,0.0],[0.0,0.0,1.0],psi)#body y
    x1= rotQuat([0.0,1.0,0.0,0.0],[0.0,0.0,1.0],psi)#intermediate
    q2= [cos(th),(yb[2:4].*sin(th))...]
    xb= rotQuat(x1,yb[2:4],th)#body x
    q3= [cos(phi),(xb[2:4].*sin(phi))...]

    #Define one quaternion for the rotation.
    qout= multQuat(multQuat(q3,q2),q1)
    magQ= norm(qout) #Make sure magnitude is 1
    if magQ==0.0
        qout = qout
    else
        qout = qout./magQ
    end
    #extract axis and angle from qout
    angle = acos(qout[1])*2
    magQ=norm(qout[2:4])
    if magQ==0.0
        axis = qout[2:4]
    else
        axis = qout[2:4]./magQ
    end

    return axis, angle
end

"""

    finite_jacobian(s, params)

    Calculates the jacobian matrix of partial derivatives for the 6-DOF model
    represented by s and params.

**Inputs**
- s = x, y, z, phi, theta, psi, u, v, w, p, q, r (same order as State)
- params = control, massproperties, reference, aeromodel, propmodel, inertialmodel, atmmodel

**Outputs**
- jacobian: a jacobian matrix representing the system.

"""

function finite_jacobian(state,params)

    #How much to step (perturb) each state value. The closer this is to zero, the more theoretically accurate.
    per=0.000001

    #Allocate ds. This makes them the same size and type. This value will be overwritten
    ds = zeros(length(state))
    ds[1:3]=state[7:9]
    #required for use of sixdof!() Not implemented.
    time = 10.0
    #Also must be reset each time sixdof!() is used.
    centerS=state
    centerP=params
    sixdof!(ds, centerS, centerP, time)#sixdof!() modifies ds so that it is now the defailt derivatives
    #Result: values in ds give the derivatives about the center point.
    centerds=ds


    #initialize the Jacobian matrix with zeros.
    jacobian=zeros(length(state),length(state))

    #Loop to perturb each element of state and build the jacobian matrix
    for i in 1:length(state)
        #Must be reset each time sixdof!() is used to ensure that we are measuring around the same point each time.
        time = 0.0
        centerS=zeros(length(state))
        centerS=state
        centerP=params
        ds = zeros(length(state))
        ds[1:3]=state[7:9]

        #Perturb s at correct index
        centerS[i]=centerS[i]+per
        #Run sixdof!() for perturbed value.
        sixdof!(ds, centerS, centerP, time)

        #Calculate the derivitives relative to centerds
        deltaS = ds - centerds

        #Store the new vector in its respective column.
        jacobian[:,i]=deltaS./per

        centerS[i]=centerS[i]-per#This may fix a weird bug
    end

    return jacobian
end

function getjacobian(s,params)
    ds = zeros(length(s))
    function centered_sixdof!(ds, s)#used exclusivley for functionality with ForwardDiff.jl
        t = 0.0
        sixdof!(ds,s,params,t)
    end
    jacobian = ForwardDiff.jacobian(centered_sixdof!, ds, s)
    return jacobian
end

"""

    geteig(s, params)

    Calculates the eigenvalues and eigenvectors for the 6-DOF model
    represented by s and params.

**Inputs**
- s = x, y, z, phi, theta, psi, u, v, w, p, q, r (same order as State)
- params = control, massproperties, reference, aeromodel, propmodel, inertialmodel, atmmodel

**Outputs**
- eigVal: Eigenvalues representing the system.
- eigVect: Eigenvectors of the system.

"""
function geteig(s, params)
    #get the jacobian matrix
    jacobian=finite_jacobian(s,params)

    eigVal, eigVect = eigen(jacobian)

    return eigVal, eigVect#eigVal(i) corresponds to eigVect(:,i)
end


"""

    animvtk(jacobian, model)

    Makes an animation of stability modes based on the jacobian matrix and a model aircraft.

**Inputs**
- jacobian: jacobian matrix for the system.
- model: Vector of arbitrary length. Each element contains info about an aerodynamic surface.
- an element of model = Wing(position, direction,  cr, ct, tilt, twist) (dir has enough info for length, sweep and dihedral)

**Outputs**
- Writes a file:

"""
function animvtk(jacobian, model, planeName)#Makes a VTK animation
    #TODO
    targetamp=pi/6#30 deg max

    #TODO: determine these two in a more robust way
    resolution = 100.0
    simEnd = 0.5

    #We need these for the animation.
    eigVal, eigVect = eigen(jacobian)
    numEig=length(eigVal)
    surfs=length(model)#Number of surfaces in the model.

    #Define necessary points for zero displacement and zero angular displacement.
    #They will be quaternions so that they can be moved easily.
    numPts=4*surfs+1
    quatPoints= zeros(4,numPts)              #4 points for each surface plus one for COM
    for i in 1:surfs                         #define all points, one surface at a time
        tilt = (pi/180)*model[i].tilt        #convert to radians
        twist = (pi/180)*model[i].twist
        quatPoints[2:4,4*i-2]= model[i].pos + [cos(tilt) 0.0 sin(tilt)].*model[i].cr/2               #Root Front
        quatPoints[2:4,4*i-1]= model[i].pos + model[i].dir + [cos(tilt+twist) 0.0 sin(tilt+twist)].*model[i].ct/2  #Tip Front
        quatPoints[2:4,  4*i]= model[i].pos - [cos(tilt) 0.0 sin(tilt)].*model[i].cr/2                             #Root Rear
        quatPoints[2:4,4*i+1]= model[i].pos + model[i].dir - [cos(tilt+twist) 0.0 sin(tilt+twist)].*model[i].cr/2  #Tip rear
    end
    rotQuatPoints=zeros(4,numPts)           #allocate
    rotQuatPoints[:,1:numPts] = quatPoints  #initialize

    #DEFINE POLYGONS
    polys= [MeshCell(PolyData.Polys(), (4*i-2):(4*i)) for i = 1:surfs]
    polys2= [MeshCell(PolyData.Polys(), (4*i-1):(4*i+1)) for i = 1:surfs]
    append!(polys,polys2)

    #CLEAR EXISTING Anim FOLDER
    rm("Anim/$planeName",recursive=true,force=true)
    mkdir("Anim/$planeName")
    mkdir("Anim/$planeName/vtk")
    numPer=0             #Counter for periodic modes
    numSub=0             #Counter for stable modes
    numUn=0              #Counter for unstable modes
    numPerSub=0          #Counter for stable periodic modes
    numPerUn=0           #Counter for unstable periodic modes

    for mode in 1:numEig #Loop for all animations

        x = eigVect[:,mode]  #get eigenvector
        lam = eigVal[mode]   #get eigenvalue
        if imag(lam)<0.0     #Modes with negative imaginary parts are duplicates
            @goto animEnd
        end

        #Label to identify what type of mode it is. Initialized as " "
        type = " "

        #initialize
        tau=1.0/real(lam)             #time constant
        TT=2.0*pi/imag(lam)           #oscilation period
        #Determine simulation length and amplitude
        if tau==NaN || tau==Inf   #Means pure periodic mode
            if TT==NaN || TT==Inf #Means Stationary mode
                simEnd=1.0
                amp=0.0
                type="Stationary"
                @goto animEnd     #Stationary modes aren't interesting
            else
                simEnd=TT          #Simulate for 1 period (if periodic)
                amp=targetamp/norm(x[4:6])
                type="Periodic"
                numPer=numPer+1
                filename="Anim/$planeName"*"/$type"*"_$numPer"
            end
        elseif TT==NaN  || TT==Inf           #Means pure exponential
            if tau<0.0                       #tau is positive (stable mode)
                simEnd=-3.0*tau              #simulate for 3 time constants
                amp=targetamp/norm(x[4:6])
                type="Subsidence"
                numSub=numSub+1
                filename="Anim/$planeName"*"/$type"*"_$numSub"
            else                            #tau is positive (unstable mode)
                simEnd=tau                  #simulate 1 time constant
                #Magnitude at end of simulation is biggest, so use it as a baseline.
                amp=targetamp/(norm(x[4:6])*real(exp(lam*simEnd)))
                type="Unstable"
                numUn=numUn+1
                filename="Anim/$planeName"*"/$type"*"_$numUn"
            end
        else                      #means exponential and periodic
            if tau<0.0            #Means decay
                simEnd=-5.0*tau   #simulate for 5 time constants
                amp=targetamp/(norm(x[4:6]))
                type="Periodic_Subsidence"
                numPerSub=numPerSub+1
                filename="Anim/$planeName"*"/$type"*"_$numPerSub"
            else                     #tau is positive (unstable mode)
                simEnd=5.0*tau       #simulate 5 time constants
                #Magnitude at end of simulation is biggest, so use it as a baseline.
                amp=2*targetamp/(norm(x[4:6])*exp(real(lam*simEnd)))
                type="Periodic_Unstable"
                numPerUn=numPerUn+1
                filename="Anim/$planeName"*"/$type"*"_$numPerUn"
            end
        end
        step=simEnd/resolution
        if amp==NaN || amp==Inf#Catch NaNs so they dont become a problem.
            amp=0.0
        end

        #ANIMATE EACH MODE
        animPVD = paraview_collection(filename*".pvd")#WARNING: Probably only works on mac and linux due to structure of "filename."

        #ADD VTK FILES TO .pvd COLLECTION. ONE VTK FILE AT A TIME
        for t in 0.0:step:simEnd

            #define state at any time, look at only real part
            xt = amp*real(x.*exp(lam*t))

            #get axis and angle from euler angles.
            rotAxis,rotAngle =eulerToQuat(xt[4:6]...)

            #rotate graphics
            for j in 1:numPts
                rotQuatPoints[:,j] = rotQuat(quatPoints[:,j],rotAxis,rotAngle)
            end

            pts = rotQuatPoints[2:4, : ]#points array (non-quaternion) just for animation
            vtkFrame = vtk_grid("Anim/$planeName"*"/vtk/mode_"*"$planeName"*"_$mode"*"__t_$t"*".vtp",pts,polys)#one vtk file
            vtkFrame["time"] = t
            vtkFrame["λreal"] = real(lam)
            vtkFrame["λimag"] = imag(lam)
            vtkFrame["τ"] = tau
            vtkFrame["T"] = TT
            vtkFrame["A"] = amp
            vtkFrame["x_x"] = real(x[1])
            vtkFrame["x_y"] = real(x[2])
            vtkFrame["x_z"] = real(x[3])
            vtkFrame["x_Φ"] = real(x[4])
            vtkFrame["x_Θ"] = real(x[5])
            vtkFrame["x_Ψ"] = real(x[6])
            vtkFrame["x_Vx"] = real(x[7])
            vtkFrame["x_Vy"] = real(x[8])
            vtkFrame["x_Vz"] = real(x[9])
            vtkFrame["x_p"] = real(x[10])
            vtkFrame["x_q"] = real(x[11])
            vtkFrame["x_r"] = real(x[12])

            animPVD[t]= vtkFrame        #paraview animation
        end
        vtk_save(animPVD)               #save the paraview animation

        #ANIMATE EACH PERIODIC MODE AGAIN, SHOWING ONLY PERIODIC PARTS
        if type=="Periodic_Unstable" || type=="Periodic_Subsidence"
            animPVD2 = paraview_collection(filename*"_Periodic"*".pvd")#WARNING: Probably only works on mac and linux
            simEnd=2*pi          #Simulate for 1 period
            step=simEnd/resolution
            amp=targetamp/norm(x[4:6])

            #ADD VTK FILES TO .pvd COLLECTION. ONE VTK FILE AT A TIME
            for t in 0.0:step:simEnd
                xt = amp*real(x).*sin(t) #define at any time, look at only real part of periodic response.

                #get axis and angle from euler angles.
                rotAxis,rotAngle =eulerToQuat(xt[4:6]...)

                #rotate graphics
                for j in 1:numPts
                    rotQuatPoints[:,j] = rotQuat(quatPoints[:,j],rotAxis,rotAngle)
                end

                pts = rotQuatPoints[2:4, : ]#points array (non-quaternion) just for animation
                vtkFrame = vtk_grid("Anim/$planeName"*"/vtk/mode_"*"$planeName"*"_$mode"*"_Periodic__t_$t"*".vtp",pts,polys)#one vtk file
                vtkFrame["time"] = t
                realLam = real(lam) #TODO Allocate these earlier to avoid doing it every loop.
                imagLam = imag(lam)
                vtkFrame["λ"] = "$realLam"*" + $imagLam"*"i"
                vtkFrame["τ"] = tau
                vtkFrame["T"] = TT
                vtkFrame["A"] = amp
                animPVD2[t*TT/(2*pi)]= vtkFrame#paraview animation
            end
            vtk_save(animPVD2)#save the paraview animation

            @label animEnd #stationary modes will skip down to here
        end

    end
end

function trimstate()
end

end # module
