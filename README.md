This package is the cubic spline interpolation for one-dimensional data. 
It is inspired by the following developments.

http://www.cs.cornell.edu/courses/cs4210/2015fa/CVLBook/new_page_1.htm


The derivative of the cubic spline is included.

"""
CSpline_coefs(x,y,EndOrder,muL,muR)

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



"""
function CSplinefDiff(x,y,EndOrder,muL,muR)
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
