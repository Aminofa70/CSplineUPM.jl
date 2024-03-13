# program utils

"""
CSpline_coefs(x,y,EndOrder1,muL,muR)
this function returns the coefficients of a cubic spline
x:: the "horizontal axis" or "abscissa."
y:: the "vertical axis" or "ordinate."
typof(x)::Vector{Float64}
typeof(x)::Vector{Float64}
EndOrder:: order of derivatives at the left and right
EndOrder:: 1 or 2; (1 is for the first order);(2 is for the second order)
typeof(EndOrder)::Int64
muL:: the value of derivative at the left side
muR:: the value of derivative at the right side
typeof(muL)::Float64
typeof(muR)::Float64
Remark: CSpline_coefs(x,y) namely without EndOrder1,muL,muR is for not-a-knot boundary condition
For example:
```
CSpline_coefs(x,y,1,0.,0.)-> clamped boundary condition
CSpline_coefs(x,y,2,0.,0.)-> natural boundary condition
CSpline_coefs(x,y)-> not-a-knot boundary condition

````
##### Example #####
````
julia> x = [0., 1., 2., 3., 4.];
julia> y = exp.(x);
julia> CC = CSpline_coefs(x,y,1,0.,0.)
julia> CC.a
4-element Vector{Float64}:
1.0
2.718281828459045
7.38905609893065
20.085536923187668
julia> CC.b
4-element Vector{Float64}:
0.0
3.941566877103152
3.400900788379342
34.556595253565355
julia> CC.c
4-element Vector{Float64}:
1.213278608273984
2.728288268829172
-3.2689543575529854
34.424648822739
julia> CC.d
4-element Vector{Float64}:
0.5050032201850616
-1.9990808754607179
12.564534393430662
-34.46863096634778
````
"""


struct SplineCoefficients
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

function CSpline_coefs(x::Vector{Float64}, y::Vector{Float64},
    EndOrder::Union{Int64, Nothing}=nothing, 
    muL::Union{Float64, Nothing}=nothing,
    muR::Union{Float64, Nothing}=nothing) :: SplineCoefficients
    n = length(x)
    ny= length(y)
    (n != ny) && error("Leng of y (=$ny) no equal to length of x (=$n)")
    (n < 3)   && error("Length of data (=$n) must be > 3")
    for i=2:n
        (x[i] <= x[i-1]) && error("x values must be monotonically increasing")
    end

    Dx = diff(x)
    yp = diff(y) ./ Dx
    T = spzeros(n-2, n-2)
    r = zeros(n-2)
    s = zeros(n)
    for i = 2:n-3
        T[i, i] = 2 * (Dx[i] + Dx[i+1])
        T[i, i-1] = Dx[i+1]
        T[i, i+1] = Dx[i]
        r[i] = 3 * (Dx[i+1] * yp[i] + Dx[i] * yp[i+1])
    end
    ########################################################
    if EndOrder == 2
        T[1, 1] = 2 * Dx[1] + 1.5 * Dx[2]
        T[1, 2] = Dx[1]
        r[1] = 1.5 * Dx[2] * yp[1] + 3 * Dx[1] * yp[2] + Dx[1] * Dx[2] * muL / 4

        T[n - 2, n - 2] = 1.5 * Dx[n - 2] + 2 * Dx[n - 1]
        T[n - 2, n - 3] = Dx[n - 1]
        r[n - 2] = 3 * Dx[n - 1] * yp[n - 2] + 1.5 * Dx[n - 2] * yp[n - 1] - Dx[n - 2] * Dx[n - 1] * muR

        stilde = T \ r

        s1 = (3 * yp[1] - stilde[1] - muL * Dx[1] / 2) / 2
        sn = (3 * yp[end] - stilde[end] + muR * Dx[end] / 2) / 2

        @inbounds s = [s1; stilde; sn]

    elseif EndOrder == 1

        T[1, 1] = 2 * (Dx[1] + Dx[2])
        T[1, 2] = Dx[1]
        r[1] = 3 * (Dx[2] * yp[1] + Dx[1] * yp[2]) - Dx[2] * muL

        T[n-2, n-2] = 2 * (Dx[n-2] + Dx[n-1])
        T[n-2, n-3] = Dx[n-1]
        r[n-2] = 3 * (Dx[n-1] * yp[n-2] + Dx[n-2] * yp[n-1]) - Dx[n-2] * muR

        @inbounds s = [muL; T \ r[1:n-2]; muR]

    elseif EndOrder === nothing && muL === nothing && muR === nothing
        q = Dx[1]^2 / Dx[2]
        T[1, 1] = 2 * Dx[1] + Dx[2] + q
        T[1, 2] = Dx[1] + q
        r[1] = Dx[2] * yp[1] + Dx[1] * yp[2] + 2 * yp[2] * (q + Dx[1])
        q = Dx[n-1]^2 / Dx[n-2]
        T[n-2, n-2] = 2 * Dx[n-1] + Dx[n-2] + q
        T[n-2, n-3] = Dx[n-1] + q
        r[n-2] = Dx[n-1] * yp[n-2] + Dx[n-2] * yp[n-1] + 2 * yp[n-2] * (Dx[n-1] + q)
        stilde = T \ r
        s1 = -stilde[1] + 2 * yp[1]
        s1 += (Dx[1] / Dx[2])^2 * (stilde[1] + stilde[2] - 2 * yp[2])
        sn = -stilde[n-2] + 2 * yp[n-1]
        sn += (Dx[n-1] / Dx[n-2])^2 * (stilde[n-3] + stilde[n-2] - 2 * yp[n-2])
        @inbounds s = [s1; stilde; sn]
    else
        error("Invalid boundary condition specified: $EndOrder")
    end
    
    a = y[1:n-1]
    b = s[1:n-1]
    Dx = diff(x)
    Dy = diff(y)
    yp = Dy ./ Dx
    c = zeros(n-1)
    d = zeros(n-1)
    @inbounds c .= (3*yp - 2*s[1:n-1]-s[2:n]) ./ (Dx )
    @inbounds d .= (s[2:n] + s[1:n-1] - 2 .* yp) ./ (Dx .* Dx)
    
    return SplineCoefficients(a, b, c, d)

end # end function cspline_coefs
##########################################################################

"""

function linspace(xstart,xend,n)
generates a linearly spaced vector
n:: is number of points
##### Example #####
````

julia> x = linspace(0.,10., 4)
4-element Vector{Float64}:
  0.0
  3.333333333333333
  6.666666666666666
 10.0

```

"""
function linspace(xstart::Float64,xend::Float64,n :: Int64)

    vector = collect(LinRange(xstart,xend,n))

    return vector::Vector{Float64}
    
end # end function linspace

##########################################################################
"""
function CSplinef(x,y,EndOrder,muL,muR)
builds the function for cubic spline
    ##### Example #####
    ````
    
    julia> xdata = [0., 1., 2.5, 3.6, 5., 7., 8.1, 10.]
    julia> ydata = sin.(xdata)
    julia> f1 = CSplinef(xdata,ydata)  # not-a-knot BC
    julia> f2 = CSplinef(xdata,ydat,1,0.,0.) #clamped BC
    julia> f3 = CSplinef(xdata,ydata,,2,0.,0.) # natural BC
    julia> f1(3)
    0.14653252176233017
    julia> f2(3)
    0.1287259900493635
    julia> f3(3)
    0.14302859118434916
    
    ```
"""

function CSplinef(x::Vector{Float64}, y::Vector{Float64},
    EndOrder::Union{Int64, Nothing}=nothing, 
    muL::Union{Float64, Nothing}=nothing,
    muR::Union{Float64, Nothing}=nothing)
    
    CC = CSpline_coefs(x, y,EndOrder, muL,muR)
    a = CC.a
    b = CC.b
    c = CC.c
    d = CC.d
    n = length(x)

    function f(xquery)
        if xquery < x[1] || xquery > x[end]
            throw(DomainError(xquery, "The query value is out of range"))
        end

        i = searchsortedfirst(x, xquery)
        i = min(max(i, 2), n) - 1
        dx = xquery - x[i]
        return a[i] + b[i] * dx + c[i] * dx^2 + d[i] * dx^3
    end
    
    return f::Function

end # end function CSplinef
##########################################################################
"""
function CSplinef(x,y,EndOrder,muL,muR)
    builds the derivative function for cubic spline
        ##### Example #####
    ````
        julia> xdata = [0., 1., 2.5, 3.6, 5., 7., 8.1, 10.]
        julia> ydata = sin.(xdata)
        julia> f = CSplinefDiff(xdata,ydata,2,0.,0.)
        julia> f(7.5)
        0.3369639804838759

    ```
"""



function CSplinefDiff(x::Vector{Float64}, y::Vector{Float64},
    EndOrder::Union{Int64, Nothing}=nothing, 
    muL::Union{Float64, Nothing}=nothing,
    muR::Union{Float64, Nothing}=nothing)
    
    CC = CSpline_coefs(x, y,EndOrder, muL,muR)
    a = CC.a
    b = CC.b
    c = CC.c
    d = CC.d
    n = length(x)

    function f(xquery)
        if xquery < x[1] || xquery > x[end]
            throw(DomainError(xquery, "The query value is out of range"))
        end

        i = searchsortedfirst(x, xquery)
        i = min(max(i, 2), n) - 1
        dx = xquery - x[i]
        return  b[i]  + 2 *c[i] * dx + 3*d[i] * dx^2
    end
    
    return f::Function

end # end function CSplinefDiff
##########################################################################

# end program uitls

