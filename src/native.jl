module NativeISEA
using StaticArrays: SVector, @SVector
using Rotations: RotX, RotZ
using GeometryOps.UnitSpherical: UnitSphereFromGeographic, GeographicFromUnitSphere, UnitSphericalPoint, slerp
const USP = UnitSphericalPoint
using LinearAlgebra: cross, dot,norm, normalize

export ISEA, ISEA20, InvISEA20, ISEA10, InvISEA10

const NORMAL_OFFSET = -0.7946544722917661229555309283275940420265905883092648017549557750084386449717329
const ITRIANGLE_SIDE_SQ = 0.90450849718747371205114670859140952943007729495144071553386215567
const ITRIANGLE_HEIGHT_SIDE = 1.044436448670983612734708092275980120462245664528624497238884055281095681373414
const MAX_TRIANGLE_HEIGHT = Float64(sind(BigFloat(60)))
const TRIANGLE_AREA = 0.6283185307179586
const TRIANGLE_PRODUCT = 0.7608452130361227

function compute_vertices()
    vertices = Vector{SVector{3,BigFloat}}(undef, 12)
    latring = atand(BigFloat(0.5))

    coslat = cosd(latring)
    sinlat = sind(latring)

    vertices[1] = @SVector [BigFloat(0.0), 0, 1]
    vertices[12] = @SVector [BigFloat(0.0), 0, -1]

    for i = 0:4
        northlon = BigFloat(- 72 * i)
        southlon = northlon - 36
        vertices[2+i] = @SVector [sind(northlon) * coslat, cosd(northlon) * coslat, sinlat]
        vertices[7+i] = @SVector [sind(southlon) * coslat, cosd(southlon) * coslat, -sinlat]
    end

    # orientation  = (BigFloat(31.7174744114611), BigFloat(78.8))
    # r1 = RotX(deg2rad(orientation[1]))
    # r2 = RotZ(deg2rad(orientation[2]))
    # rot = r2 * r1

    #rvertices = (rot,).*vertices
    rvertices = vertices
    return UnitSphericalPoint.(rvertices)
end

function A(a, b, c)
    x = triple_product(a, b, c)
    y = (1 + dot(a, b) + dot(b, c) + dot(c, a))
    if abs(y) < 1e-10
        return 2 * atan(x / y), false
    end
    # @show (1 + dot(a, b) + dot(b, c) + dot(c, a))
    2 * atan(x / y),true 
end
triple_product(a, b, c) = dot(a, cross(b, c))

struct ISEATriangle{T}
    a::UnitSphericalPoint{T}
    b::UnitSphericalPoint{T}
    c::UnitSphericalPoint{T}
    kind::Symbol
    #And store normal and offset of the corresponding plane
    # normal::UnitSphericalPoint{T}
    # offset::T
end
function ISEATriangle(a, b, c, kind, T)
    # normal = -cross(a-b,c-b)
    # normal = normal/norm(normal)
    # offset = -dot(normal,a)
    return ISEATriangle(USP{T}(a),USP{T}(b),USP{T}(c),kind)
end
function Base.show(io::IO, ::MIME"text/plain", tri::ISEATriangle)
    corners = GeographicFromUnitSphere().((tri.a, tri.b, tri.c))
    println(io, "$(tri.kind) triangle with corner coordinates $corners")
end
_base(t::ISEATriangle) = t.c-t.b
_height(t::ISEATriangle) = t.a-(t.b+t.c)/2


#Make a type for the neighboring triangles
struct ISEANeighbor
    i::Int #Index of the neighboring triangle
    type::Int # Type of the connection, 0:left/right, same direction, 1:left/right, opposite direction, 2:up/down
end
function makeneighbors()
    #Make neighbors for north triangles
    northneighbors = map(1:5) do i
        ISEANeighbor(mod1(i-1,5),0),ISEANeighbor(mod1(i+1,5),0),ISEANeighbor(i+15,2)
    end
    southneighbors = map(6:10) do i
        ISEANeighbor(mod1(i+1,5)+5,0),ISEANeighbor(mod1(i-1,5)+5,0),ISEANeighbor(21-i,2)
    end
    northtipneighbors = map(11:15) do i
        ISEANeighbor(mod1(i-1,5)+15,1),ISEANeighbor(i+5,1),ISEANeighbor(21-i,2)
    end
    southtipneighbors = map(16:20) do i
        ISEANeighbor(i-5,1),ISEANeighbor(mod1(i+1,5)+10,1),ISEANeighbor(i-15,2)
    end
    return vcat(northneighbors,southneighbors,northtipneighbors,southtipneighbors)
end

"""
    ISEA{T}

A structure to hold the ISEA triangles and their neighbors.
"""
struct ISEA{T}
    triangles::Vector{ISEATriangle{T}}
    neighbors::Vector{NTuple{3,ISEANeighbor}}
end
function ISEA(T=Float64)
    vertices = compute_vertices()
    northindices = [(1,2,3),(1,3,4),(1,4,5),(1,5,6),(1,6,2)]
    southindices = [(12,11,10),(12,10,9),(12,9,8),(12,8,7),(12,7,11)]
    northtriangles = [ISEATriangle(vertices[i], vertices[j], vertices[k], :north, T) for (i, j, k) in northindices]
    southtriangles = [ISEATriangle(vertices[i], vertices[j], vertices[k], :south, T) for (i, j, k) in southindices]
    # Generate the triangles for the equator, we make sure that the triangles
    equatorindices_northtip = [(2,11,7),(3,7,8),(4,8,9),(5,9,10),(6,10,11)]
    equatorindices_southtip = [(7,3,2),(8,4,3),(9,5,4),(10,6,5),(11,2,6)]
    equatorindices = vcat(equatorindices_northtip,equatorindices_southtip)
    equatortriangles = [ISEATriangle(vertices[i], vertices[j], vertices[k], :middle, T) for (i, j, k) in equatorindices]

    triangles = vcat(northtriangles,southtriangles,equatortriangles)
    neighbors = makeneighbors()
    ISEA(triangles,neighbors)
end


struct ISEA20{T} <: Function
    isea::ISEA{T}
end
struct InvISEA20{T} <: Function 
    isea::ISEA{T}
end
Base.inv(isea::ISEA20) = InvISEA20(isea.isea)
Base.inv(isea::InvISEA20) = ISEA20(isea.isea)
ISEA20(args...;kwargs...) = ISEA20(ISEA(args...;kwargs...))
InvISEA20(args...;kwargs...) = InvISEA20(ISEA(args...;kwargs...))

(isea::ISEA20)(latlon::NTuple{2}) = isea(UnitSphericalPoint(latlon))
function _transform_isea(isea::ISEA,p::UnitSphericalPoint)
    grid = isea
    r, stable = transform_point(p, grid.triangles[1])
    d, _ = triangle_distance(r)
    stable && iszero(d) && return (1, r...)
    current_min = d, (1, r...)
    for i in 2:20
        r,stable = transform_point(p,grid.triangles[i])
        d, _ = triangle_distance(r)
        stable && iszero(d) && return (i, r...)
        if d < first(current_min)
            current_min = d, (i, r...)
        end
    end
    #Nothing was found, so lets use the triangle with the smallest distance
    return last(current_min)
end

(isea::ISEA20)(p::UnitSphericalPoint) = _transform_isea(isea.isea,p)

# function transform_point(p::UnitSphericalPoint,tri::ISEATriangle)
#     t = -tri.offset/dot(p,tri.normal)
#     t < 0 && return (Inf,Inf)
#     pointonsurface = p*t
#     nbase=_base(tri)*ITRIANGLE_SIDE_SQ
#     nheight = _height(tri) * ITRIANGLE_HEIGHT_SIDE
#     x1 = dot(pointonsurface-tri.b,nbase)
#     x2 = dot(pointonsurface-tri.b,nheight)
#     x1,x2
# end

@inline function transform_bary(v0,v1,v2,v)

    p1 = TRIANGLE_PRODUCT*v - triple_product(v,v1,v2)*v0
    all(iszero, p1) || (p1 = normalize(p1))
    #dot(p1,v1) < 0 && (p1 = -p1)
    h = sqrt((1 - dot(v0, v)) / (1 - dot(v0, p1)))
    a_part, stable = A(v0,v1,p1)
    β2 = h * a_part / TRIANGLE_AREA
    β0 = 1 - h
    β1 = h-β2
    β0, β1, β2, stable
end
function bary_to_xy(β0, _, β2)
    y = 0.5 * sqrt(3) * β0
    x = 0.5 * β0 + β2
    x, y
end

function transform_point(v, t)
    β0, β1, β2, stable = transform_bary(t.a,t.b,t.c,v)
    bary_to_xy(β0, β1, β2), stable
end

#intriangle((x1, x2)) = (0 <= x1 <= 1) && (0 <= x2 <= MAX_TRIANGLE_HEIGHT * (1 - 2 * abs(x1 - 0.5)))
intriangle((x1, x2)) = iszero(first(triangle_distance((x1, x2))))
function intriangle(p,tri::ISEATriangle) 
    (x,y),stable = transform_point(UnitSphericalPoint(p),tri)
    stable && intriangle((x,y))
end

function triangle_distance((x1, x2))
    downdist = max(-x2, zero(x2))
    x1offset = x1 - 0.5
    leftdist = max(x2 - MAX_TRIANGLE_HEIGHT * (1 + 2 * (x1offset)), zero(x2))
    rightdist = max(x2 - MAX_TRIANGLE_HEIGHT * (1 - 2 * (x1offset)), zero(x2))
    findmax((leftdist, rightdist, downdist))
end


# function itransform_point(x1,x2,tri::ISEATriangle)
#     base=_base(tri)
#     height=_height(tri)
#     pointonsurface = tri.b + x1 * base + x2 * height / MAX_TRIANGLE_HEIGHT
#     p = pointonsurface/norm(pointonsurface)
#     return GeographicFromUnitSphere()(p)
# end
function itransform_point(x, y, t)
    β0 = 2 * y / sqrt(3)
    β2 = x - 0.5 * β0
    h = 1 - β0
    iszero(h) && return t.a
    q = if β2 != 0.0
        a = β2 / h * TRIANGLE_AREA
        S = sin(a)
        C = 1 - cos(a)
        #@show a,S,C
        f = S * TRIANGLE_PRODUCT + C * (dot(t.a, t.b) * dot(t.b, t.c) - dot(t.c, t.a))
        g = C * sqrt(1 - dot(t.b, t.c)^2) * (1 + dot(t.a, t.b))
        #@show f,g
        2 / acos(dot(t.b, t.c)) * atan(g / f)
    else
        0.0
    end
    #@show q
    p = slerp(t.b, t.c, q)
    T = acos(1 + h^2 * (dot(t.a, p) - 1)) / acos(dot(t.a, p))
    slerp(t.a, p, T)
end
(isea::InvISEA20)((n, x1, x2)) = itransform_point(x1, x2, isea.isea.triangles[n])


northpairs = [(i,i+15) for i in 1:5]
southpairs = [(i,21-i) for i in 11:15]
diamond2tri = vcat(northpairs,southpairs)
tri2diamond = [(1,1),(2,1),(3,1),(4,1),(5,1),
               (10,2),(9,2),(8,2),(7,2),(6,2),
               (6,1),(7,1),(8,1),(9,1),(10,1),
               (1,2),(2,2),(3,2),(4,2),(5,2)]
struct ISEA10{T}<:Function
    isea::ISEA{T}
end
ISEA10(args...;kwargs...) = ISEA10(ISEA(args...;kwargs...))
function _twenty_to_ten(n,x1,x2)
    i,j = tri2diamond[n]
    if j == 2
        #We are in the south diamond, so we need to reverse the x2 coordinate
        x1 = 1-x1
        x2 = -x2
    end
    return i,x1,x2
end
function (isea::ISEA10)(p::UnitSphericalPoint)
    _twenty_to_ten(_transform_isea(isea.isea,p)...)
end
(isea::ISEA10)(latlon::NTuple{2}) = isea(UnitSphericalPoint(latlon))

struct InvISEA10{T} <: Function
    isea::ISEA{T}
end
Base.inv(isea::ISEA10) = InvISEA10(isea.isea)
Base.inv(isea::InvISEA10) = ISEA10(isea.isea)
function _ten_to_twenty(i,x1,x2)
    i1,i2 = diamond2tri[i]
    i,x1,x2 = x2 < 0 ? (i2,1-x1,-x2) : (i1,x2)
    i,x1,x2
end
function (isea::InvISEA10)((i, x1, x2)) 
    i,x1,x2 = _ten_to_twenty(i,x1,x2)
    return itransform_point(x1,x2,isea.isea.triangles[i])
end

struct ISEA5{T} <: Function
    isea::ISEA{T}
end
ISEA5(args...;kwargs...) = ISEA5(ISEA(args...;kwargs...))
rect2diamond(i) = (i,i+5)
diamond2rect(i) = mod1(i,5),i÷5+1
function (isea::ISEA5)(p::UnitSphericalPoint)
    i,x1,x2 = _twenty_to_ten(_transform_isea(isea.isea,p)...)
    i,lr = diamond2rect(i)
    if lr == 2
        x1 = 0.5-x1
        x2 = -sqrt(3)/2-x2
    end
    i,x1,x2
end
(isea::ISEA5)(latlon::NTuple{2}) = isea(UnitSphericalPoint(latlon))




end