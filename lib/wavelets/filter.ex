defmodule Wavelets.Filter do
  @moduledoc """
  Wavelet filter coefficients and related utilities.
  """

  @type t :: %__MODULE__{
          name: String.t(),
          family: atom(),
          vanishing_moments: integer() | nil,
          decomposition_low_pass: list(float()),
          decomposition_high_pass: list(float()),
          reconstruction_low_pass: list(float()),
          reconstruction_high_pass: list(float()),
          support_width: pos_integer()
        }

  defstruct [
    :name,
    :family,
    :vanishing_moments,
    :decomposition_low_pass,
    :decomposition_high_pass,
    :reconstruction_low_pass,
    :reconstruction_high_pass,
    :support_width
  ]
end
