from resilience_measurement import run_model_multi

seq = [5, (2, 5), (10, 14), 12, (7, 8)]   # buses are ints; road links are (u, v)

_ = run_model_multi(
    sequence=seq,
    result_folder="Experiment/2025-11-11_13-00-00",
    message="Multi-crew demo",
    Scenario="SENS1",          # uses your file-copy logic
    plot_control=True,
    focus=False,
    power_crews=2,             # strictly power
    road_crews=3,              # strictly road
    service_time_power=20.0,   # matches your original bus repair time
    service_time_road=10.0     # matches your original link repair time
)