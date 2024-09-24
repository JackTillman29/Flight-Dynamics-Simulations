function outVect = ReflectVector(inVect,uNorm)
    outVect = 2 * uNorm * dot(inVect,uNorm) - inVect;
end