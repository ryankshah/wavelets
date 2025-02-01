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
  def transform_1d(signal, filter, levels, _precision \\ :double) do
    do_transform_1d(signal, filter, levels, [], [])
  end

  defp do_transform_1d(signal, _filter, 0, approx_acc, detail_acc) do
    {Enum.reverse([signal | approx_acc]), Enum.reverse(detail_acc)}
  end

  defp do_transform_1d(signal, filter, level, approx_acc, detail_acc) do
    {approx, details} = undecimated_decomposition(signal, filter)
    next_filter = upsample_filter(filter)

    do_transform_1d(
      approx,
      next_filter,
      level - 1,
      [approx | approx_acc],
      [details | detail_acc]
    )
  end

  defp undecimated_decomposition(signal, filter) do
    # Extend signal by periodic padding to maintain length
    n = length(signal)
    padded = extend_periodic(signal, n)

    # Apply filters without downsampling
    approx = convolve_periodic(padded, filter.decomposition_low_pass, n)
    details = convolve_periodic(padded, filter.decomposition_high_pass, n)

    {approx, details}
  end

  defp extend_periodic(signal, n) do
    signal ++ Enum.take(signal, n)
  end

  defp convolve_periodic(signal, filter, output_length) do
    0..(output_length - 1)
    |> Enum.map(fn i ->
      filter
      |> Enum.with_index()
      |> Enum.map(fn {f, k} ->
        # Use periodic indexing
        idx = rem(i + k, output_length)
        Enum.at(signal, idx) * f
      end)
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
    Enum.flat_map(coefficients, &[&1, 0.0])
  end
end
