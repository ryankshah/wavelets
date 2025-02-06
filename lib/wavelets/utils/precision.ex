defmodule Wavelets.Precision do
  @moduledoc """
  Handles precision-specific computations for wavelet transforms
  """

  @doc """
  Rounds a number to the appropriate precision
  """
  def round_to_precision(value, :single) do
    <<f::float-32>> = <<value |> Float.round(7)::float-32>>
    f
  end

  def round_to_precision(value, :double), do: value

  @doc """
  Applies precision to a list of numbers
  """
  def apply_precision(list, precision) when is_list(list) do
    Enum.map(list, &round_to_precision(&1, precision))
  end

  @doc """
  Gets machine epsilon for the given precision
  """
  # 2^-23
  def epsilon(:single), do: 1.1920929e-7
  # 2^-52
  def epsilon(:double), do: 2.220446049250313e-16

  @doc """
  Checks if two values are equal within precision
  """
  def approx_equal?(a, b, precision) do
    abs(a - b) < epsilon(precision)
  end
end
