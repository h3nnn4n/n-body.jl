include("types.jl")
include("config.jl")

function main()
    particles = Array(particle, npar)

    for i in 1:npar
        particles[i] = particle(rand(), [rand() * lbox, rand() * lbox, rand() * lbox], [rand(), rand(), rand()], [0.0, 0.0, 0.0])
    end

    t  = 0.0f0
    a  = [0.0f0, 0.0f0, 0.0f0] :: Array{Float32, 1}
    rj = [0.0f0, 0.0f0, 0.0f0] :: Array{Float32, 1}

    while t < tmax
        for i in 1:npar
            for ic in 1:3
                print("$(particles[i].position[ic]) ")
                a[ic]  = 0
                rj[ic] = 0
            end
            print("\n")

            for j in 1:npar
                if i != j
                    for ic in 1:3
                        rj[ic] = particles[j].position[ic]
                    end

                    dist  = 0.0f0
                    for ic in 1:3
                        dist += particles[i].position[ic]^2 + rj[ic]^2
                    end
                    dist = dist^(0.5f0)
                    dist = (dist*dist + epss * epss)^(0.5)
                    for ic in 1:3
                        a[ic] += -GC*particles[j].mass/(dist^3)*( particles[i].position[ic] - rj[ic] )
                    end
                end
            end

            for ic in 1:3
                particles[i].position[ic]     = particles[i].position[ic] + particles[i].velocity[ic]*tstp + 0.5*particles[i].acceleration[ic]*tstp*tstp
                particles[i].velocity[ic]     = particles[i].velocity[ic] + 0.5*( particles[i].acceleration[ic] + a[ic] )*tstp
                particles[i].acceleration[ic] = a[ic]
            end
        end

        t += tstp
    end
end

main()
