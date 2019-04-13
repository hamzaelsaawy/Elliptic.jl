#
# Tests for Landen Sequences
#

@testset "Landen" begin
    using Elliptic.Landen

    areseqequal(x, y; args...) = all(isapprox.(x, y; args...))

    @testset "basics" begin
        @test LandenSeq(1, descending=true) isa NonConvergedLandenSeq
        @test LandenSeq(0, descending=false) isa NonConvergedLandenSeq
        @test LandenSeq(-eps(0.0)) isa NonConvergedLandenSeq
        @test LandenSeq(1 + eps()) isa NonConvergedLandenSeq
        @test LandenSeq(NaN) isa NonConvergedLandenSeq
        @test LandenSeq(NaN, descending=false) isa NonConvergedLandenSeq
        @test LandenSeq((0.5, ), descending=false) isa LandenSeq{0, Base.Bottom}
        @test LandenSeq("0.5", descending=false) isa LandenSeq{0, Base.Bottom}

        l0 = LandenSeq(0)
        l1 = LandenSeq(1)

        @test l0.ks == l1.k′s == (0.0, )
        @test l1.ks == l0.k′s == (1.0, )
    end

    @testset "$(f)" for f in ["table1", "table2"]
        ks = first( readdlm(joinpath("data/papers/ow", "$f.csv"), ',', header=true) )[:, 1]
        k′s = Landen._k′.(ks)
        # data given to 9 decimal places
        rtol_data = 1e-9

        @testset "$T" for T in [Float32, Float64, BigFloat]
            rtol_T = eps(T)
            rtol = max(rtol_T, rtol_data)

            k = T(ks[1])
            k′ = T(k′s[1])
            # paper aims for 16 digits of accuaracy
            lc = LandenSeq(k, descending=true, ktol=1e-16)
            lc′ = LandenSeq(k′, descending=false, ktol=1e-16)

            # test _k′
            @test areseqequal(hypot.(lc.ks, lc.k′s), 1, atol=rtol_T())
            @test areseqequal(hypot.(lc′.ks, lc′.k′s), 1, atol=rtol_T())
            # descending seq is correct
            @test areseqequal(lc.ks, ks, rtol=1e-9, atol=rtol_T())
            @test areseqequal(lc.k′s, k′s, rtol=1e-9, atol=rtol_T())
            # ascending seq works
            @test areseqequal(lc′.k′s, ks, rtol=1e-9, atol=rtol_T())
            @test areseqequal(lc′.ks, k′s, rtol=1e-9, atol=rtol_T())

            N = length(ks)-1
            @test LandenSeq(k, N=N-1, descending=true, ktol=1e-16) isa NonConvergedLandenSeq
        end
    end



    @testset "a&b table 17.1" begin
        # read as bigfloat since Float64 to BigFloat introduces errors on order of 1e-17
        t171, _ = readdlm(joinpath("data/ab", "table_17_1.csv"), ',', BigFloat, header=true)
        # data given to 15th decimal place
        atol_data = 1e-15

        ms = t171[:, 1]
        Ks = t171[:, 2]
        N = size(t171, 1)

        @testset "$T" for T in [Float32, Float64, BigFloat]
            atol_T = eps(T)
            rtol_T = eps(T) * 3

            @testset for i in 1:N
                ii = N - i + 1
                KK = Ks[i]
                KK′ = Ks[ii]
                k = sqrt(T(ms[i]))
                k′ = sqrt(T(ms[ii]))

                l = LandenSeq(k, k′, ktol=100*atol_T)
                @test areseqequal(hypot.(l.ks, l.k′s), 1, atol=atol_T)

                Kl, K′l = @inferred K(l)
                @test Kl isa T
                @test K′l isa T

                @test Kl ≈ KK atol=atol_data rtol=rtol_T
                @test K′l ≈ KK′ atol=atol_data rtol=rtol_T
            end
        end
    end # table 17.1

    # errors are 6× larger than above
    # likely due to using sin(α) vs computing with α directly...
    @testset "a&b table 17.2" begin
        # read as bigfloat since Float64 to BigFloat introduces errors on order of 1e-17
        t172, _ = readdlm(joinpath("data/ab", "table_17_2.csv"), ',', BigFloat, header=true)
        # data given to 15th decimal place
        atol_data = 1e-15

        αs = t172[:, 1]
        Ks = t172[:, 2]
        N = size(t172, 1)

        @testset "$T" for T in [Float32, Float64, BigFloat]
            atol_T = eps(T)
            rtol_T = eps(T) * 3

            @testset for i in 1:N
                ii = length(Ks) - i + 1
                # avoid inaccuracies from sind() by doing first as BigFloat
                k = T(sind(αs[i]))
                k′ = T(sind(αs[ii]))
                KK = Ks[i]
                KK′ = Ks[ii]

                l = LandenSeq(k, k′, ktol=100*atol_T)
                @test areseqequal(hypot.(l.ks, l.k′s), 1, atol=atol_T)

                Kl, K′l = @inferred K(l)
                @test Kl isa T
                @test K′l isa T

                @test Kl ≈ KK atol=6*atol_data rtol=rtol_T
                @test K′l ≈ KK′ atol=6*atol_data rtol=rtol_T
            end
        end
    end # table 17.2
end # end landen tests
