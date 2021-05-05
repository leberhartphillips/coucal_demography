library(RMark)
data(dipper, package = "RMark")

dipper.proc <- 
  RMark::process.data(data = dipper,
                      model = "CJS",
                      groups = "sex")
dipper.ddl <- 
  RMark::make.design.data(dipper.proc)

dipper_CJS_global <- 
  RMark::mark(model = "CJS", groups = "sex", 
              model.name = "dipper_CJS_test", filename = "dipper_CJS_test",
              data = dipper.proc, 
              ddl = dipper.ddl,
              model.parameters = list(Phi = list(formula = ~sex * Time),
                                      p = list(formula = ~Time)))
RMark::export.MARK(x = dipper.proc, project.name = "dipper_test_export", 
                   model = dipper_CJS_global, replace = TRUE)