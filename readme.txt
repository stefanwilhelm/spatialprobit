model.frame(formula, ...)      # generic method; returns a data.frame with variables needed to use 'formula' and ... arguments
model.frame.lm()               # specific to linear models
model.frame.glm()              # specific to glm models

model.frame.lm <- function (formula, ...) 
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1L]] <- as.name("lm")
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env)) 
            env <- parent.frame()
        eval(fcall, env, parent.frame())
    }
    else formula$model
}

model.frame.glm <- function (formula, ...) 
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1L]] <- as.name("glm")
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env)) 
            env <- parent.frame()
        eval(fcall, env)
    }
    else formula$model
}


####################################################################################################

model.matrix()            # creates a design (or model) matrix.

model.matrix.lm <- function (object, ...) 
{
    if (n_match <- match("x", names(object), 0L)) 
        object[[n_match]]
    else {
        data <- model.frame(object, xlev = object$xlevels, ...)
        NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    }
}




####################################################################################################

model.response(data, type = "any") # returns the response of a model frame
model.weights()                    # returns the weights of a model frame
model.offset()                     # Include an Offset in a Model Formula; An offset is a term to be added to a linear predictor, such as in a generalised linear model, with known coefficient 1 rather than an estimated coefficient. 
model.extract()

