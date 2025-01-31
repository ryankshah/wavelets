defmodule Wavelets.IDWT do
  @moduledoc """
  Implementation of the Inverse Discrete Wavelet Transform
  """

  alias Wavelets.Filter

  @doc """
  Performs 1D inverse discrete wavelet transform
  """
  @spec inverse_1d(list(number), list(number), Filter.t(), Wavelets.precision()) ::
          list(number)
  def inverse_1d(approximation, details, filter, _precision \\ :double) do
    # Scale coefficients back
    scale = :math.sqrt(2)
    scaled_approx = Enum.map(approximation, &(&1 * scale))
    scaled_details = Enum.map(details, &(&1 * scale))

    # Upsample and apply reconstruction filters
    upsampled_approx =
      upsample_convolve(scaled_approx, filter.reconstruction_low_pass)

    upsampled_details =
      upsample_convolve(scaled_details, filter.reconstruction_high_pass)

    # Combine and trim padding
    Enum.zip(upsampled_approx, upsampled_details)
    |> Enum.map(fn {a, d} -> a + d end)
    |> trim_signal(filter.support_width)
  end

  @doc """
  Performs 2D inverse discrete wavelet transform
  """
  @spec inverse_2d(
          list(list(number)),
          {list(list(number)), list(list(number)), list(list(number))},
          Filter.t(),
          Wavelets.precision()
        ) :: list(list(number))
  def inverse_2d(
        approximation,
        {h_details, v_details, d_details},
        filter,
        precision \\ :double
      ) do
    # Inverse transform on columns
    rows_low =
      transpose(approximation)
      |> Enum.zip(transpose(v_details))
      |> Enum.map(fn {a, d} -> inverse_1d(a, d, filter, precision) end)
      |> transpose()

    rows_high =
      transpose(h_details)
      |> Enum.zip(transpose(d_details))
      |> Enum.map(fn {h, d} -> inverse_1d(h, d, filter, precision) end)
      |> transpose()

    # Inverse transform on rows
    Enum.zip(rows_low, rows_high)
    |> Enum.map(fn {l, h} -> inverse_1d(l, h, filter, precision) end)
  end

  # Helper functions
  defp upsample_convolve(signal, filter) do
    # Upsample
    upsampled = signal |> Enum.flat_map(&[&1, 0.0])

    # Add padding for convolution
    padding_size = length(filter) - 1

    padded =
      List.duplicate(0.0, padding_size) ++
        upsampled ++ List.duplicate(0.0, padding_size)

    # Convolve
    0..(length(upsampled) + length(filter) - 2)
    |> Enum.map(fn i ->
      Enum.zip(
        Enum.slice(padded, i..(i + length(filter) - 1)),
        filter
      )
      |> Enum.map(fn {s, f} -> s * f end)
      |> Enum.sum()
    end)
  end

  defp trim_signal(signal, filter_length) do
    padding = div(filter_length - 1, 2)
    len = length(signal)
    Enum.slice(signal, padding..(len - padding - 1))
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end
end
