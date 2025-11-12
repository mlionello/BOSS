function resU = get_residual(U, infeat)
    bU = infeat\U;
    resU = U-infeat*bU;
end
