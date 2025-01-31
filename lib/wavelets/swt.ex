defmodule Wavelets.SWT do
  @moduledoc """
  Implementation of the Stationary Wavelet Transform (Undecimated Wavelet Transform)
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D stationary wavelet transform
  """
  @spec transform_1d(
          list(number),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) ::
          {list(list(number)), list(list(number))}
  def transform_1d(signal, filter, levels, precision \\ :double) do
    do_transform_1d(signal, filter, levels, [], [])
  end

  defp do_transform_1d(signal, _filter, 0, approx_acc, detail_acc) do
    {Enum.reverse([signal | approx_acc]), Enum.reverse(detail_acc)}
  end

  defp do_transform_1d(signal, filter, level, approx_acc, detail_acc) do
    # Apply filters without downsampling
    {approx, details} = undecimated_decomposition(signal, filter)

    # Upsample filter for next level
    next_filter = upsample_filter(filter)

    do_transform_1d(approx, next_filter, level - 1, [approx | approx_acc], [
      details | detail_acc
    ])
  end

  # Helper functions
  defp undecimated_decomposition(signal, filter) do
    extended = extend_signal_swt(signal, length(filter.decomposition_low_pass))

    approx = convolve(extended, filter.decomposition_low_pass)
    details = convolve(extended, filter.decomposition_high_pass)

    {approx, details}
  end

  defp extend_signal_swt(signal, filter_length) do
    padding_size = filter_length - 1
    padding = Enum.take(signal, padding_size)
    padding ++ signal ++ Enum.reverse(padding)
  end

  defp convolve(signal, filter) do
    filter_length = length(filter)

    0..(length(signal) - filter_length)
    |> Enum.map(fn i ->
      Enum.zip(Enum.slice(signal, i..(i + filter_length - 1)), filter)
      |> Enum.map(fn {s, f} -> s * f end)
      |> Enum.sum()
    end)
  end

  defp upsample_filter(%Filter{} = filter) do
    %{
      filter
      | decomposition_low_pass:
          upsample_coefficients(filter.decomposition_low_pass),
        decomposition_high_pass:
          upsample_coefficients(filter.decomposition_high_pass),
        reconstruction_low_pass:
          upsample_coefficients(filter.reconstruction_low_pass),
        reconstruction_high_pass:
          upsample_coefficients(filter.reconstruction_high_pass),
        support_width: filter.support_width * 2
    }
  end

  defp upsample_coefficients(coefficients) do
    Enum.flat_map(coefficients, &[&1, 0])
  end
end
