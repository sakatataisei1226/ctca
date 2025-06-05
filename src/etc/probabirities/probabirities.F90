module m_phase_probabilities
    use m_probability_base, only: t_PhaseProbability
    use m_probability_maxwell, only: t_MaxwellWithVelocityProbability, new_MVProbability
    use m_probability_grid, only: t_GridPhaseProbability, new_GridPhaseProbability
    implicit none

    private
    public t_PhaseProbability

    public t_MaxwellWithVelocityProbability
    public new_MVProbability

    public t_GridPhaseProbability
    public new_GridPhaseProbability

end module
