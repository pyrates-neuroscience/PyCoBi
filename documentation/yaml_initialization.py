from pycobi import ODESystem

a = ODESystem.from_yaml("model_templates.neural_mass_models.qif.qif",
                        auto_dir="", init_cont=False)

