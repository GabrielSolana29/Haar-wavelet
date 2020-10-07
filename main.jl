using Plots,DataFrames,CSV

## Functions unrelated to wavelets

function samplingData(size)
     x = 0
     y = zeros(size)
     for i =1:size
         x = i
         y[i] = 2*(x^2) - 3*x + 1
     end
     display(plot(1:size,y[1:size]))
     #gui(plot(1:size,y[1:size]))
     return y
end

### if wav is 1 it would return the matrix with one and - one
function one_matrix_filter(size,size_coef,wav)
    rows = Int(size/2)
    mat = zeros(rows,Int(size))
    cont = 1
    for i = 1:rows
        for j = 1: Int(size_coef/rows)
            if wav ==1 && j%2 ==0
                mat[i,cont] = -1
                cont = cont + 1
            else
                mat[i,cont] = 1
                cont = cont + 1
            end
        end
    end
    return mat
end

function energy_magnitude_aproximation(a,d,currentLevel,noLevels)
    x,y = size(a)
    sum = 0
    for i = 1:x
        sum = sum + (a[i,currentLevel+1])^2
    end
    return sum
end

function energy_magnitude_aproximation_original(sampled_signal)
    x = size(sampled_signal)
    sum_a = 0
    for i = 1:x[1]
        sum_a = sum_a + (sampled_signal[i])^2
    end
    return sum_a
end

function energy_magnitude_details(a,d,currentLevel,noLevels)
    x,y = size(d)
    sum = 0
    #for l = 2:currentLevel+1
        for i = 1:x
            sum = sum + (d[i,currentLevel+1])^2
        end
    #end
    return sum
end

function energy_levels(a,d,noLevels,complete)
    energy_vect = zeros(noLevels,2)
    #complete = energy_magnitude_aproximation(a,d,0,noLevels)
    for i =1:noLevels
        energy_vect[i,1] = (energy_magnitude_aproximation(a,d,i,noLevels) * 100)/complete
        energy_vect[i,2] = (energy_magnitude_details(a,d,i,noLevels) * 100)/complete
    end
    display(plot(1:noLevels,energy_vect[:,1],title = "Aproximation energy",xaxis = "Level of decomposition", yaxis = "Energy Percentage",label="aproximation",legend =:bottomleft)) #legend=false
    display(plot!(1:noLevels,energy_vect[:,2],title = "Details energy",xaxis = "Level of decomposition", yaxis = "Energy Percentage",label="details"))
    return energy_vect
end

## Function realted to wavelets

function sup_val(level,support)
    x = sqrt(1/support)
    return x
end

function approximation_coeff(no_samples,level,h_coef,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    ## Create the ones and zeros matrix and multiply with the filter coefficients to obtain the filter matrix
    mat = one_matrix_filter(size_prev_level,size_prev_level,0)
    ## If analysis == 1 matrix is not transpose, synthesis is performed otherwise
    filter_mat = h_coef * mat
    ## Extract the value of the aproximation coefficients from previous level
    approximation = a[1:size_prev_level,level]
    ## Multiply filter with previous approximation
    filter_coeff = filter_mat * approximation
    ## Concatenate the matrix to be the same size as the array a
    filter_mat_complete = vcat(filter_coeff,zeros(no_samples-size_coeff))
    return filter_mat_complete
end

function wavelet_coeff(no_samples,level,h_coef,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    ## Create the ones and zeros matrix and multiply with the filter coefficients to obtain the filter matrix
    mat = one_matrix_filter(size_prev_level,size_prev_level,1)
    filter_mat = h_coef * mat
    ## Extract the value of the aproximation coefficients from previous level
    approximation = a[1:size_prev_level,level]
    ## Multiply filter with previous approximation
    filter_coeff = filter_mat*approximation
    ## Concatenate the matrix to be the same size as the array a
    filter_mat_complete = vcat(filter_coeff,zeros(no_samples-size_coeff))

    return filter_mat_complete
end

function synthesis_approximation_coeff(no_samples,level,h_coef,a2,d)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    ## Create the ones and zeros matrix and multiply with the filter coefficients to obtain the filter matrix
    mat_aprox = one_matrix_filter(size_prev_level,size_prev_level,0)
    mat_detail = one_matrix_filter(size_prev_level,size_prev_level,1)
    filter_mat_aprox = h_coef * mat_aprox
    filter_mat_detail = h_coef * mat_detail
    ## Extract the value of the aproximation coefficients from next level
    approximation = a2[1:size_coeff,1]
    ## Multiply filter with previous approximation (the matrix is transposed)
    filter_coeff = (filter_mat_aprox' * approximation) + (filter_mat_detail' * d[1:size_coeff,1])
    x = size(filter_coeff)
    ## Concatenate the matrix to be the same size as the array a
    filter_mat_complete = vcat(filter_coeff,zeros(no_samples-x[1]))

    return filter_mat_complete
end

## Wavelet Analysis Function

function waveletAnalysis(noLevels,no_samples,sampled_signal)

    ## For the Haar transform the number of coefficients is 2. h = (1/sqrt(2)),(1/sqrt(2)) and h1 = (1/sqrt(2)), -(1/sqrt(2))
    size_coef = 2
    h_coef = 1/sqrt(2)
    ## Zeros matrices with approximation coefficients and wavelets coefficients
    a = zeros(no_samples,noLevels+1)
    d = zeros(no_samples,noLevels+1)
    # First column of the matrix is the sampled signal or level 0 of decomposition (julia starts the indexes in 1 instead of 0)
    a[:,1] = sampled_signal'
    ## Obtain coefficients for each level
    for i = 1:noLevels+1
        if i != noLevels+1
            ## Calculate the coefficients
            a[:,i+1] = approximation_coeff(no_samples,i,h_coef,a)

            d[:,i+1] = wavelet_coeff(no_samples,i,h_coef,a)
        end
    end

    return a,d
end

function waveletSynthesis(noLevels,no_samples,sampled_signal,a,d,currentLevel)
    ## For the Haar transform the number of coefficients is 2. h = (1/sqrt(2)),(1/sqrt(2)) and h1 = (1/sqrt(2)), -(1/sqrt(2))
    size_coef = 2
    h_coef = 1/sqrt(2)
    ## Zeros matrices with approximation coefficients and wavelets coefficients
    a2 = zeros(no_samples,noLevels+1)
    d2 = zeros(no_samples,noLevels+1)
    a2[1,currentLevel+1] = a
    ## Obtain coefficients for each level
    for i = currentLevel+1:-1:1
        if i != noLevels+1
            ## Calculate the coefficients
            a2[:,i] = synthesis_approximation_coeff(no_samples,i,h_coef,a2[:,i+1],d[:,i+1])
            #return a2
        end
    end

    return a2
end

function plot_aproximations(a,d,no_samples,noLevels)
    x,y = size(a)
    a_recon = zeros(x,y)
    a_recon[:,1] = a[:,1]
    d_recon = zeros(x,y)
    cont  = noLevels
    for i = 1:Int(noLevels)
        a_1 = zeros(no_samples)
        support =  Int(2^i)
        ## Obtain the value for the support if the approximated signals is reconstructed with the coeff
        val_support= sup_val(i,support)
        sp = zeros(support)
        sp .= val_support
        sp_d = zeros(support)
        for z = 1:Int(support/2)
            sp_d[z] = val_support
        end

        for z = Int(support/2)+1:support
            sp_d[z] = -val_support
        end

        new = sp * a[1,i+1]
        new_d = sp_d * d[1,i+1]
        for j = 2:Int(no_samples/(2^i))
            hlp = sp * a[j,i+1]
            hlp_d = sp_d * d[j,i+1]
            new = vcat(new,hlp)
            new_d = vcat(new_d,hlp_d)
        end
        a_recon[:,i+1] = new
        d_recon[:,i+1] = new_d
    end

    ## this needs to be displayed but with the reconstructed aproximation
    siz = 1/(6)
    l3 = @layout grid(noLevels+1, 1, heights=[siz,siz,siz,siz])
    l10 = @layout grid(6, 2, heights=[siz,siz,siz,siz,siz,siz])
    #display(plot(1:no_samples,a_recon[:,1:noLevels+1], layout = (noLevels+1, 1),xlims = (0,no_samples),ylims = (0,maximum(a[:,1]))))
    #display(plot(1:no_samples,a_recon[:,1:noLevels+1], layout = l10,xlims = (0,no_samples),ylims = (0,maximum(a[:,1]))))
    # plot details
    display(plot(1:no_samples,d_recon[:,1:noLevels+1], layout = l10,xlims = (0,no_samples),ylims = (-1000000,maximum(a[:,1]))))
    return a_recon,d_recon
end

## Function for performing compression
function compressionAprox_coef(a,d,level,no_samples,energy_threshold)
    ## Obtain the enrgy of the original signal
    ax,ay = size(a)
    maxEnergy = 0
    for i = 1:ax
         #global maxEnergy = maxEnergy + (a[i,1])^2
         maxEnergy = maxEnergy + (a[i,1])^2
    end
    # Find the energy of each coefficient
    no_samples_level = Int(no_samples / (2^level))
    newVec = zeros(no_samples_level)
    for i = 1:no_samples_level
        #global level
        #global newVec[i] = (a[i,level+1])^2
         newVec[i] = (a[i,level+1])^2
    end
    ## this method returns the order by the smallest to the largest of the vector
    index_vec = sortperm(newVec)
    ## Eliminate half of the coefficients
    #respaldo = copy(newVec)
    half_coeff =  copy(newVec)
    energy_threshold_vec = copy(newVec)

    ## Find the compression by eleiminating half of the smallest coefficients
    for i = 1:no_samples_level
        if index_vec[i] <= no_samples_level/2
            half_coeff[i] = 0
        end
    end
    ## find compression by setting a threshold
    sum_energy = 0
    cont = no_samples_level
    cont_2= 0
    sum_energy_percentage = 0
    while sum_energy_percentage < energy_threshold

        if cont == 0
            print("\n Threshold couldn't be reached!!")
            return 0
        end
        ## To avoid infinite looping
        ## iterates in reverse because it goes from the larges to the smallest coefficient
        sum_energy = sum_energy + energy_threshold_vec[index_vec[cont]]
        sum_energy_percentage = ((sum_energy * 100) / maxEnergy)
        cont_2 = cont_2 + 1

        cont = cont - 1
    end

    for i = 1:no_samples_level
        if i > cont_2
            energy_threshold_vec[index_vec[no_samples_level + 1 - i]] = 0
        end
    end

    energy_threshold_percentage = (sum(energy_threshold_vec)*100)/maxEnergy
    energy_threshold_vec = sqrt.(energy_threshold_vec)
    half_coeff_percentage = (sum(half_coeff) * 100) / maxEnergy
    half_coeff = sqrt.(half_coeff)
    print("\n The number of samples required to reach the threshold are: ", cont_2)
    print("\n The energy percentage obtained by removing half the coefficients is: ", half_coeff_percentage)
    print("\n The energy percentage obtained by using the energy threshold is: ", sum_energy_percentage)

    return half_coeff_percentage,half_coeff,energy_threshold_percentage,energy_threshold_vec, cont_2
end

function reconstruct_aproximations(comp,no_samples,level)
    a_recon = zeros(no_samples)
    a_1 = zeros(no_samples)
    support =  Int(2^level)
    ## Obtain the value for the support if the approximated signals is reconstructed with the coeff
    val_support= sup_val(level,support)
    sp = zeros(support)
    sp .= val_support
    new = sp * comp[1]
    for j = 2:Int(no_samples/(2^level))
        hlp = sp * comp[j]
        new = vcat(new,hlp)
    end
    return new
end


## Execute the main program
no_samples = 1024
noLevels = 10
### Obtain the samples from the function and plot it
sampled_signal = samplingData(no_samples)
a,d = @time waveletAnalysis(noLevels,no_samples,sampled_signal)
current_level = noLevels
a_reconstruction = @time waveletSynthesis(noLevels,no_samples,sampled_signal,a[1,current_level+1],d,current_level)
### Obtain the remaining energy on the approximation
complete = energy_magnitude_aproximation_original(sampled_signal)
energy_vect = energy_levels(a,d,noLevels,complete)
### plot reconstruction of aproximation
a_recon,d_recon  = plot_aproximations(a[:,:],d,no_samples,noLevels)
### Perform compression with the approximation coefficients
energy_threshold = 99;
level = 6
# Returns the perecentage of energy that removing half the coefficients generate, and a vector with the coefficients that remained
# It returns the percentage of energy that removing x nuimber of coefficients resulted, and the vector with the coefficients, as well as the number of coefficients needed
half_c_percentage,half_c,threshold_percntg,threshold_vec,no_samples_for_compression = compressionAprox_coef(a,d,level,no_samples,energy_threshold)
# reconstruct the compressed signals
plot(1:Int(no_samples/2^level),a[1:Int(no_samples/2^level),level])
comp = half_c

cs = zeros(Int(no_samples/2^level),1)
cs .= comp
CSV.write("half_c.csv",  DataFrame(cs), writeheader=false)
plot(1:(no_samples/2^level),comp)

reconstruction = reconstruct_aproximations(comp,no_samples,level)

### Add the compressed with details from previous levels
for i = 1:level
    global reconstruction
    reconstruction = reconstruction + d_recon[:,level + 1-i]
end


cs = zeros(Int(no_samples),1)
cs .= reconstruction
CSV.write("reconstruction_half_c.csv",  DataFrame(cs), writeheader=false)
plot(1:no_samples,reconstruction)
# now reconstruct the compressed signal that used the threshold of energy

comp = threshold_vec
cs = zeros(Int(no_samples/2^level),1)
cs .= comp
CSV.write("threshold_energy.csv",  DataFrame(cs), writeheader=false)
plot(1:(no_samples/2^level),comp)


reconstruction = reconstruct_aproximations(comp,no_samples,level)
### Add the compressed with details from previous levels
for i = 1:level
    global reconstruction
    reconstruction = reconstruction + d_recon[:,level + 1-i]
end

cs = zeros(no_samples,1)
cs .= reconstruction
CSV.write("reconstruction_threshold_energy.csv",  DataFrame(cs), writeheader=false)
plot(1:no_samples,reconstruction)
#plot(1:(no_samples/2^level),threshold_vec)

## Save in csv file the results of the reconstruction
CSV.write("Scaling_coeff.csv",  DataFrame(a), writeheader=false)
CSV.write("Wavelet_coeff.csv",  DataFrame(d), writeheader=false)
CSV.write("energy.csv",  DataFrame(energy_vect), writeheader=false)
CSV.write("aprox_reconstruction.csv",  DataFrame(a_recon), writeheader=false)
CSV.write("details_reconstruction.csv",  DataFrame(d_recon), writeheader=false)


#=
## Test with signal from exercise
sampled_signal = [4 6 10 12 8 6 5 5]
no_samples = 8
noLevels = 3
a,d = @time waveletAnalysis(noLevels,no_samples,sampled_signal)
current_level = noLevels
a_reconstruction = @time waveletSynthesis(noLevels,no_samples,sampled_signal,a[1,current_level+1],d,current_level)
### Obtain the remaining energy on the approximation
complete = energy_magnitude_aproximation_original(sampled_signal)
energy_vect, = energy_levels(a,d,noLevels,sampled_signal)
### plot reconstruction of aproximation
a_recon,d_recon = plot_aproximations(a,d,no_samples,noLevels)

energy_threshold = 90;
level = 1
# Returns the perecentage of energy that removing half the coefficients generate, and a vector with the coefficients that remained
# It returns the percentage of energy that removing x nuimber of coefficients resulted, and the vector with the coefficients, as well as the number of coefficients needed
half_coeff_percentage,half_coeff,energy_threshold_percentage,energy_threshold_vec,no_samples_for_compression = compressionAprox_coef(a,d,level,no_samples,energy_threshold)
=#
