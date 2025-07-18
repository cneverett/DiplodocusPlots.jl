""" 
    SIUnitsTime(t_code)

Returns the time in SI units seconds from the time in code units of seconds * σTc 
"""
function SIUnitsTime(t_code::T) where T <: AbstractFloat
    # Convert time code to SI units (seconds)
    σT = getfield(DC,Symbol("σT"))
    c = getfield(DC,Symbol("c"))  # σT is the time unit in seconds

    t_SI::T = t_code / (σT*c)  # Convert time code to seconds
    return t_SI
end
function SIUnitsTime()
    return raw"$[\text{s}]$"
end

""" 
    CodeUnitsTime(t_code)

Returns the time in code units from the time in code units
"""
function CodeUnitsTime(t_code::T) where T <: AbstractFloat
    return t_code
end
function CodeUnitsTime()
    return raw"$[\text{s} \times \sigma_{T}c]$"
end
