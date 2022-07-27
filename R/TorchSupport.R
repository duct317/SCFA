#' @import torch
SCFA_dataset <- dataset(

    name = "SCFA_dataset",

    initialize = function(input) {
        self$data <- self$prepare_scDHA_data(input)
    },

    .getbatch = function(index) {
        x <- self$data[index, ]
        list(x)
    },

    .length = function() {
        self$data$size()[[1]]
    },

    prepare_scDHA_data = function(input) {
        input <- as.matrix(input)
        torch_tensor(input)
    }
)

SCFA_AE <- nn_module(
    "SCFA_AE",
    initialize = function(original_dim, im_dim) {
        self$fc1 <- nn_linear(original_dim, im_dim)
        self$fc2 <- nn_linear(im_dim, original_dim)

        nn_init_zeros_(self$fc1$bias)
        nn_init_zeros_(self$fc2$bias)
        nn_init_xavier_uniform_manual(self$fc1$weight)
        nn_init_xavier_uniform_manual(self$fc2$weight)
    },
    forward = function(x) {
        x %>%
        self$fc1() %>%
        self$fc2()
    }
)

nn_init_xavier_uniform_manual <- function(tensor, gain = 1) {
    fans <- nn_init_calculate_fan_in_and_fan_out(tensor)
    fan_in <- fans[[1]]
    fan_out <- fans[[2]]
    std <- gain * sqrt(2.0 / (fan_in + fan_out))
    a <- sqrt(3.0) * std # Calculate uniform bounds from standard deviation
    nn_init_no_grad_uniform(tensor, -a, a)
}

nn_init_calculate_fan_in_and_fan_out <- function(tensor) {
    dimensions <- tensor$dim()
    num_input_fmaps <- tensor$size(2)
    num_output_fmaps <- tensor$size(1)
    receptive_field_size <- 1

    # if (dimensions > 2)
    #   receptive_field_size <- tensor[1,1,..]$numel()

    fan_in <- num_input_fmaps * receptive_field_size
    fan_out <- num_output_fmaps * receptive_field_size

    list(fan_in, fan_out)
}

nn_init_no_grad_uniform <- function(tensor, a, b) {
    with_no_grad({
        out <- tensor$uniform_(a, b)
    })
    out
}

state <- function(self) {
    attr(self, "state")
}

`state<-` <- function(self, value) {
    attr(self, "state") <- value
    self
}

optim_adamw <- optimizer(
    "optim_adamw",
    initialize = function(params, lr=1e-3, betas=c(0.9, 0.999), eps=1e-8,
                            weight_decay=1e-2) {
        if (lr < 0)
            value_error("Invalid learning rate: {lr}")

        if (eps < 0)
            value_error("Invalid eps: {eps}")

        if (betas[[1]] < 0 || betas[[1]] > 1)
            value_error("Invalid beta parameter at index 1")

        if (betas[[2]] < 0 || betas[[2]] > 1)
            value_error("Invalid beta parameter at index 2")

        if (weight_decay < 0)
            value_error("Invalid weight decay value: {weight_decay}")

        defaults <- list(lr=lr, betas=betas, eps = eps, weight_decay = weight_decay)

        super$initialize(params, defaults)
    },

    step = function(closure = NULL) {
        loop_fun <- function(group, param, g, p) {
            grad <- param$grad

            # state initialization
            if (length(state(param)) == 0) {
                state(param) <- list()
                state(param)[["step"]] <- 0
                state(param)[["exp_avg"]] <- torch_zeros_like(param, memory_format=torch_preserve_format())
                state(param)[["exp_avg_sq"]] <- torch_zeros_like(param, memory_format=torch_preserve_format())
            }

            # Perform stepweight decay
            param$mul_(1 - group$lr*group$weight_decay)

            exp_avg <- state(param)[["exp_avg"]]
            exp_avg_sq <- state(param)[["exp_avg_sq"]]

            beta1 <- group$betas[[1]]
            beta2 <- group$betas[[2]]

            state(param)[["step"]] <- state(param)[["step"]] + 1
            bias_correction1 <- 1 - beta1 ^ state(param)[['step']]
            bias_correction2 <- 1 - beta2 ^ state(param)[['step']]

            # Decay the first and second moment running average coefficient
            exp_avg$mul_(beta1)$add_(grad, alpha=1 - beta1)
            exp_avg_sq$mul_(beta2)$addcmul_(grad, grad, value=1 - beta2)

            denom <- (exp_avg_sq$sqrt() / sqrt(bias_correction2))$add_(group$eps)

            step_size <- group$lr / bias_correction1

            param$addcdiv_(exp_avg, denom, value=-step_size)
        }
        private$step_helper(closure, loop_fun)
    }
)

check_grad_nan <- function(parameters) {
    if (is(parameters, "torch_tensor"))
        parameters <- list(parameters)
    parameters <- Filter(function(x) !is_undefined_tensor(x$grad), parameters)
    for (p in parameters) {
        if(p$grad$sum()$isnan()$item()) return(TRUE)
    }
    FALSE
}
