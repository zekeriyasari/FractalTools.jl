# This file includes tests for IFS

@testset "IFSTestset" begin 
    #------------------ Test transformation construction ---------------------------------------- #

    # Check the fields
    for name in [:A, :b] 
        @test hasfield(Transformation, name)
    end
    
    # Check Transformation construction 
    A, b = rand(2,2), rand(2)
    w = Transformation(A, b)
    @test w.A == A 
    @test w.b == b 

    # Check dimension and contraction factor
    @test dimension(w) == size(A, 2)
    @test contfactor(w) ≈ norm(A) 

    # Dimension mismatch in w(x) = A*x + b
    @test_throws DimensionMismatch Transformation(rand(2,2), rand(3))  

    # --------------------  Test IfS construction ----------------------------------------------- #  

    # Check fields 
    for name in [:ws, :probs] 
        @test hasfield(IFS, name)
    end

    # Check IFS construction 
    ws = [Transformation(rand(2,2), rand(2)) for i in 1 : 4]
    probs = [0.1, 0.2, 0.3, 0.4]
    ifs = IFS(ws, probs)
    @test ifs.ws == ws 
    @test ifs.probs == probs 

    # Check dimension and contraction factors 
    @test dimension(ifs) == 2
    @test contfactor(ifs) == maximum(contfactor.(ifs.ws))

    # If not specified, the probs of transformations are equal.
    ifs = IFS(ws) 
    @test ifs.ws == ws 
    @test ifs.probs ≈ 1 / length(ws) *  ones(length(ws))
   
    # Dimension of the transformations must match
    @test_throws DimensionMismatch IFS([
        Transformation(rand(3,3), rand(3)),     # 3-dimensional transformation 
        Transformation(rand(3,3), rand(3)),     # 3-dimensional transformation 
        Transformation(rand(2,2), rand(2))      # 2-dimensional transformation 
        ]
    )

    # --------------------  Test attractor construction ----------------------------------------------- #  
    
    # Attractor of IFS for line attractor 
    ifs = IFS([
        Transformation(fill(0.5, 1, 1), fill(0., 1)), 
        Transformation(fill(0.5, 1, 1), fill(1 / 2, 1))
    ]) 
    atr = attractor(ifs, [[0.25]])
    @test atr.alg == DetAlg()
    @test atr.ifs == ifs 
    @test atr.initset == [[0.25]]
    @test atr.parallel == false
    @test atr.numiter == 10

    initset = [rand(1) for i in 1 : 5]
    atr = attractor(ifs, initset, alg=DetAlg(), numiter=5, parallel=true)
    @test atr.alg == DetAlg() 
    @test atr.numiter = 5 
    @test atr.parallel = false 
end
