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
        sum = sum + a[i,currentLevel+1]^2
    end
    return sum
end

function energy_magnitude_details(a,d,currentLevel,noLevels)
    x,y = size(a)
    sum = 0
    for l = 2:currentLevel+1
        for i = 1:x
            sum = sum + d[i,l]^2
        end
    end
    return sum
end

function energy_levels(a,d,noLevels)
    energy_vect = zeros(noLevels,2)
    complete = energy_magnitude_aproximation(a,d,0,noLevels)
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

## Test with signal from exercise
sampled_signal = [4 6 10 12 8 6 5 5]
no_samples = 8
noLevels = 3
a,d = @time waveletAnalysis(noLevels,no_samples,sampled_signal)
current_level = noLevels
a_reconstruction = @time waveletSynthesis(noLevels,no_samples,sampled_signal,a[1,current_level+1],d,current_level)
### Obtain the remaining energy on the approximation
energy_vect = energy_levels(a,d,noLevels)
### plot reconstruction of aproximation
a_recon,d_recon = plot_aproximations(a,d,no_samples,noLevels)

## Execute the main program
no_samples = 1024
noLevels = 10
### Obtain the samples from the function and plot it
sampled_signal = samplingData(no_samples)
a,d = @time waveletAnalysis(noLevels,no_samples,sampled_signal)
current_level = noLevels
a_reconstruction = @time waveletSynthesis(noLevels,no_samples,sampled_signal,a[1,current_level+1],d,current_level)
### Obtain the remaining energy on the approximation
energy_vect = energy_levels(a,d,noLevels)
### plot reconstruction of aproximation
a_recon,d_recon  = plot_aproximations(a[:,:],d,no_samples,noLevels)

## Save in csv file the results of the reconstruction
CSV.write("Scaling_coeff.csv",  DataFrame(a), writeheader=false)
CSV.write("Wavelet_coeff.csv",  DataFrame(d), writeheader=false)
CSV.write("energy.csv",  DataFrame(energy_vect), writeheader=false)
CSV.write("aprox_reconstruction.csv",  DataFrame(a_recon), writeheader=false)
CSV.write("details_reconstruction.csv",  DataFrame(d_recon), writeheader=false)
