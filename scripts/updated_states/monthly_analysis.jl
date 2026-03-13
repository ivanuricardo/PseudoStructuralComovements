using DrWatson
@quickactivate :PseudoStructuralComovements
using CSV, Tables, SparseArrays
using XLSX
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

# Without using Michigan or Wisconsin
# Not sure if I should still keep this series
rawdata = XLSX.readdata(datadir("./state_indexes/reguib_northcentral.xlsx"), "Sheet1!A2:S459")
vecdata = Float64.(rawdata[:, 2:end])'
coincident = vecdata[1:2:end, 232:end]'

# Before transformation
ut_employment = load_series("employment")
ut_unemployment = load_series("unemployment")
ut_hours = load_series("hours")

employment = transform(ut_employment)
unemployment = transform(ut_unemployment; type = "diff")
hours = transform(ut_hours)

file = XLSX.readxlsx(datadir("updated_states/wages.xlsx"))
sheet = file["Table"]
ut_wages = Float64.(XLSX.getdata(sheet, "C7:BJ15"))'
# Rearrange to be alphabetical and align with other series
# ia, il, in, mi, mn, nd, oh, sd, wi
rearranged_wages = ut_wages[:, [3,1,2,4,5,6,7,8,9]]
monthly_wages = denton(rearranged_wages, ut_employment)
# monthly_wages = quarterly_to_monthly(rearranged_wages)
wages = transform(monthly_wages)

ser = 1:9
Plots.plot(demean_standardize(ut_employment[:, ser]))
Plots.plot(demean_standardize(unemployment[:, ser]))
Plots.plot!(demean_standardize(wages[:, ser]))
Plots.plot(ut_hours[:, 3])
Plots.plot(coincident[:, ser])

catted_data = cat(wages', unemployment', hours', employment'; dims = 3)
tendata = permutedims(catted_data, (1,3,2))
# Corresponds to the perm WI, MI, SD, MN, ND, IA, OH, IL, IN
pp = [9, 4, 8, 5, 6, 1, 7, 2, 3]
# pp = [2,3,6]
# or perm can be IA, IL, MN, OH, SD, ND, IL
# pp = [1,2,4,6,7,3,5]
rearranged_tendata = tendata[pp, :, :]
matdata = vectorize(rearranged_tendata)

n1, n2, obs = size(rearranged_tendata)
dimvals = [n1, n2]
icest = rank_selection(matdata, dimvals; iters=1000, pmax=1)
# aic selects (7,4) with one lag
# bic selects (7,3) with one lag
# hqc selects (7,3) with one lag
# this is 100 total parameters to estimate

ser = 5
Plots.plot(hcat(employment[:, ser], wages[:, ser]))

# Wages are determined by unemployment (positive), hours (negative), and employment (positive)
res = comovement_reg(matdata, dimvals, [5, 1]; iters=1000, p=1)
other_res = rrmar(matdata, dimvals, [5,1])
# Nothing is significant except the relation between Indiana and Ohio (neighbors)
# I could e.g., restrict the second dimension to be 1. Then I get two significant coefs
# however, if I restrict the second dimension to be 1, then the first dimension
# should be 2 or 3 (2 for HQC 3 for AIC, 2 for BIC)
#
#  employment unemployment hours wages
# julia> res.gamma_est
# 4×3 Matrix{Float64}:
#   1.0        0.0      0.0
#   0.0        1.0      0.0
#   0.0        0.0      1.0
#  -0.599247*  58.6196*  -0.257812
#
# julia> res.delta_est
# 7×5 Matrix{Float64}:
#   1.0        0.0      0.0        0.0         0.0
#   0.0        1.0      0.0        0.0         0.0
#   0.0        0.0      1.0        0.0         0.0
#   0.0        0.0      0.0        1.0         0.0
#   0.0        0.0      0.0        0.0         1.0
#  -0.404534  -5.5722  -1.2895*    -0.664706   -0.723722
#  -0.54262   10.9298   0.423614   0.0167974   1.00337
#
# Monthly values from Denton method
# julia> icest.ictable'
# 28×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -53179.2  -52987.4  -53101.8  1.0  1.0  1.0
#  -53358.4  -53128.9  -53265.8  2.0  1.0  1.0 **
#  -53369.6  -53109.3  -53264.6  3.0  1.0  1.0 *
#  -53365.1  -53080.8  -53250.4  4.0  1.0  1.0
#  -53365.0  -53063.6  -53243.4  5.0  1.0  1.0
#  -53360.6  -53048.9  -53234.8  6.0  1.0  1.0
#  -53535.3  -53220.3  -53408.2  7.0  1.0  1.0
#  -53349.6  -53140.7  -53265.3  1.0  2.0  1.0
#  -53621.5  -53374.9  -53522.0  2.0  2.0  1.0
#  -53741.0  -53463.6  -53629.1  3.0  2.0  1.0
#  -53858.2  -53556.8  -53736.6  4.0  2.0  1.0
#  -53953.7  -53635.2  -53825.2  5.0  2.0  1.0
#  -54038.3  -53709.5  -53905.6  6.0  2.0  1.0
#  -54070.0  -53737.8  -53936.0  7.0  2.0  1.0
#  -53361.3  -53142.1  -53272.8  1.0  3.0  1.0
#  -53639.6  -53382.8  -53536.0  2.0  3.0  1.0
#  -53815.1  -53527.4  -53699.0  3.0  3.0  1.0
#  -53950.3  -53638.7  -53824.6  4.0  3.0  1.0
#  -54074.6  -53745.9  -53942.0  5.0  3.0  1.0
#  -54159.2  -53820.1  -54022.4  6.0  3.0  1.0
#  -54221.6  -53879.1  -54083.4  7.0  3.0  1.0 BIC (100 par)
#  -53361.4  -53138.8  -53271.5  1.0  4.0  1.0
#  -53644.4  -53384.1  -53539.4  2.0  4.0  1.0
#  -53840.3  -53549.2  -53722.9  3.0  4.0  1.0
#  -53966.4  -53651.3  -53839.2  4.0  4.0  1.0
#  -54080.2  -53748.0  -53946.2  5.0  4.0  1.0
#  -54162.4  -53819.9  -54024.2  6.0  4.0  1.0
#  -54223.0  -53877.1  -54083.4  7.0  4.0  1.0 AIC HQC (101 par)

# Monthly repeated values
# julia> icest.ictable'
# 28×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -51578.3  -51386.5  -51500.9  1.0  1.0  1.0
#  -51768.2  -51538.8  -51675.6  2.0  1.0  1.0 **
#  -51777.5  -51517.2  -51672.5  3.0  1.0  1.0 *
#  -51771.5  -51487.3  -51656.8  4.0  1.0  1.0
#  -51763.2  -51461.8  -51641.6  5.0  1.0  1.0
#  -51762.5  -51450.8  -51636.7  6.0  1.0  1.0
#  -51761.0  -51445.9  -51633.9  7.0  1.0  1.0
#  -51602.8  -51393.9  -51518.5  1.0  2.0  1.0
#  -51786.1  -51539.5  -51686.6  2.0  2.0  1.0
#  -51818.9  -51541.5  -51706.9  3.0  2.0  1.0
#  -51826.9  -51525.5  -51705.3  4.0  2.0  1.0
#  -51837.5  -51519.0  -51709.0  5.0  2.0  1.0
#  -51837.8  -51509.0  -51705.2  6.0  2.0  1.0
#  -51840.5  -51508.3  -51706.4  7.0  2.0  1.0
#  -51602.7  -51383.5  -51514.2  1.0  3.0  1.0
#  -51791.4  -51534.5  -51687.7  2.0  3.0  1.0
#  -51863.8  -51576.1  -51747.7  3.0  3.0  1.0
#  -51877.2  -51565.5  -51751.4  4.0  3.0  1.0
#  -51894.1  -51565.3  -51761.5  5.0  3.0  1.0
#  -51895.8  -51556.7  -51758.9  6.0  3.0  1.0
#  -51896.8  -51554.3  -51758.6  7.0  3.0  1.0
#  -51601.8  -51379.2  -51511.9  1.0  4.0  1.0
#  -51793.3  -51533.0  -51688.3  2.0  4.0  1.0
#  -51887.0  -51595.9  -51769.5  3.0  4.0  1.0 BIC (85 par)
#  -51899.0  -51583.9  -51771.9  4.0  4.0  1.0
#  -51910.7  -51578.5  -51776.6  5.0  4.0  1.0 HQC (97 par)
#  -51911.9  -51569.4  -51773.7  6.0  4.0  1.0 AIC (100 par)
#  -51911.6  -51565.7  -51772.0  7.0  4.0  1.0

# For 9 states
#
# julia> icest.ictable'
# 36×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -53696.3     -53450.9     -53686.8     1.0  1.0  1.0
#  -54143.2     -53849.9     -54065.1     2.0  1.0  1.0
#  -54185.2     -53850.5     -54047.6     3.0  1.0  1.0
#  -54201.9     -53832.2     -54014.1     4.0  1.0  1.0
#  -54201.4     -53803.0     -53972.4     5.0  1.0  1.0
#  -54130.7     -53710.0     -53869.7     6.0  1.0  1.0
#  -54198.6     -53762.0     -53914.8     7.0  1.0  1.0
#  -54256.7     -53810.4     -53959.1     8.0  1.0  1.0
#  -54307.0     -53857.6     -54004.9     9.0  1.0  1.0
#  -54061.8     -53836.8     -54061.8     1.0  2.0  1.0
#  -54337.5     -54028.3     -54232.5     2.0  2.0  1.0
#  -54453.1     -54102.5     -54288.6     3.0  2.0  1.0
#  -54545.0     -54159.3     -54330.2     4.0  2.0  1.0
#  -54637.4     -54223.0     -54381.5     5.0  2.0  1.0
#  -54725.4     -54288.7     -54437.5     6.0  2.0  1.0
#  -54810.6     -54358.0     -54499.8     7.0  2.0  1.0
#  -54867.6     -54405.5     -54543.1     8.0  2.0  1.0
#  -54895.0     -54429.6     -54565.9     9.0  2.0  1.0
#  -54106.4     -53835.5     -54053.9     1.0  3.0  1.0
#  -54348.7     -54029.9     -54227.5     2.0  3.0  1.0
#  -54496.5     -54136.3     -54315.9     3.0  3.0  1.0
#  -54619.5     -54224.3     -54388.6     4.0  3.0  1.0
#  -54738.1     -54314.2     -54466.0     5.0  3.0  1.0
#  -54838.7     -54392.5     -54534.7     6.0  3.0  1.0
#  -54938.0     -54475.9     -54611.1     7.0  3.0  1.0
#  -55006.0     -54534.3     -54665.3     8.0  3.0  1.0
#  -55054.6     -54579.7     -54709.4     9.0  3.0  1.0
#  -54108.7     -53834.5     -54050.7     1.0  4.0  1.0
#  -54350.1     -54028.2     -54223.6     2.0  4.0  1.0
#  -54516.1     -54152.7     -54330.1     3.0  4.0  1.0
#  -54636.2     -54237.8     -54400.0     4.0  4.0  1.0
#  -54750.5     -54323.4     -54473.0     5.0  4.0  1.0
#  -54843.7     -54394.3     -54534.2     6.0  4.0  1.0
#  -54938.7     -54473.4     -54606.4     7.0  4.0  1.0
#  -55005.3     -54530.4     -54659.3     8.0  4.0  1.0
#  -55053.2     -54575.0     -54702.5     9.0  4.0  1.0

# We don't get any reduction, except if we use the EBIC and restrict the second dimension to be one.

# With simple extrapolation
# we select 5,4 with our BIC and 2,1 with EBIC
#
# julia> icest.ictable'
# 36×6 adjoint(::Matrix{Float64}) with eltype Float64:
#  -52421.9  -52176.5  -52412.5  1.0  1.0  1.0
#  -52608.0  -52314.7  -52529.9  2.0  1.0  1.0
#  -52646.7  -52312.0  -52509.1  3.0  1.0  1.0
#  -52666.3  -52296.6  -52478.5  4.0  1.0  1.0
#  -52664.6  -52266.1  -52435.6  5.0  1.0  1.0
#  -52663.4  -52242.7  -52402.4  6.0  1.0  1.0
#  -52654.3  -52217.6  -52370.4  7.0  1.0  1.0
#  -52650.3  -52204.1  -52352.7  8.0  1.0  1.0
#  -52649.3  -52199.8  -52347.1  9.0  1.0  1.0
#  -52436.3  -52175.0  -52400.0  1.0  2.0  1.0
#  -52615.4  -52306.2  -52510.4  2.0  2.0  1.0
#  -52664.1  -52313.4  -52499.6  3.0  2.0  1.0
#  -52707.5  -52321.8  -52492.7  4.0  2.0  1.0
#  -52744.9  -52330.5  -52489.0  5.0  2.0  1.0
#  -52755.1  -52318.5  -52467.2  6.0  2.0  1.0
#  -52753.0  -52300.4  -52442.2  7.0  2.0  1.0
#  -52755.1  -52293.0  -52430.6  8.0  2.0  1.0
#  -52755.3  -52289.9  -52426.2  9.0  2.0  1.0
#  -52440.3  -52169.4  -52387.8  1.0  3.0  1.0
#  -52623.2  -52304.4  -52502.1  2.0  3.0  1.0
#  -52685.2  -52325.0  -52504.6  3.0  3.0  1.0
#  -52736.5  -52341.2  -52505.6  4.0  3.0  1.0
#  -52780.7  -52356.8  -52508.6  5.0  3.0  1.0
#  -52795.0  -52348.7  -52490.9  6.0  3.0  1.0
#  -52796.8  -52334.6  -52469.8  7.0  3.0  1.0
#  -52807.8  -52336.1  -52467.2  8.0  3.0  1.0
#  -52807.9  -52332.9  -52462.6  9.0  3.0  1.0
#  -52439.8  -52165.6  -52381.9  1.0  4.0  1.0
#  -52621.2  -52299.3  -52494.7  2.0  4.0  1.0
#  -52700.5  -52337.1  -52514.5  3.0  4.0  1.0
#  -52748.5  -52350.1  -52512.3  4.0  4.0  1.0
#  -52792.9  -52365.8  -52515.5  5.0  4.0  1.0
#  -52805.3  -52355.9  -52495.8  6.0  4.0  1.0
#  -52811.1  -52345.7  -52478.7  7.0  4.0  1.0
#  -52813.9  -52339.0  -52467.8  8.0  4.0  1.0
#  -52813.3  -52335.2  -52462.7  9.0  4.0  1.0
