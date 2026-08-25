
function halfspace_vx(x, z, V::Function)
    integrand(t) = (2 / pi) * sin(t)^2 * V(x - z * tan(t))
    val, _ = quadgk(integrand, -pi/2, pi/2)
    return val
end

function halfspace_vz(x, z, V::Function)
    integrand(t) = (2 / pi) * sin(t) * cos(t) * V(x - z * tan(t))
    val, _ = quadgk(integrand, -pi/2, pi/2)
    return val
end
