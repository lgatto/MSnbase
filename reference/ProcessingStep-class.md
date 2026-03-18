# Simple processing step class

The `ProcessingStep` class is a simple object to encapsule all relevant
information of a data analysis processing step, i.e. the function name
and all arguments.

## Details

Objects of this class are mainly used to record all possible processing
steps of an
[`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)
object for later *lazy execution*.

## Objects from the Class

Objects can be created by calls of the form `new("ProcessingStep",...)`
or using the `ProcessingStep` constructor function.

## Slots

- `FUN`::

  The function name to be executed as a character string.

- `ARGS`::

  A named `list` with all arguments to the function.

## Methods and functions

- executeProcessingStep(object, ...):

  Execute the processing step `object`. Internally this calls `do.call`
  passing all arguments defined in the `ProcessingStep` `object` along
  with potential additional arguments in `...` to the function
  `object@FUN`.

## Extends

Class
`"`[`Versioned`](https://rdrr.io/pkg/Biobase/man/class.Versioned.html)`"`,
directly.

## Author

Johannes Rainer \<johannes.rainer@eurac.edu\>

## See also

[`OnDiskMSnExp`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.md)

## Examples

``` r
## Define a simple ProcessingStep
procS <- ProcessingStep("sum", list(c(1, 3, NA, 5), na.rm= TRUE))

executeProcessingStep(procS)
#> [1] 9
```
