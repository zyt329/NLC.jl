# global mu_B, k_B, N_A

function define_constants(;mu_B_self_defined = 9.27401007828e-24, k_B_self_defined = 1.380649e-23, N_A_self_defined = 6.02214076e23)
    const mu_B = mu_B_self_defined  # joules per tesla
    const k_B = k_B_self_defined  # joules per K
    const N_A = N_A_self_defined  # Avogadro's number
    println(mu_B)
end

define_constants()
